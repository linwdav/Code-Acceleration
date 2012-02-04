//Galactic Habitable Zone code by Mike Gowanlock
//Do not distribute or publish any of this code without first obtaining permission from Mike Gowanlock.
//Date: July 11, 2011


#include <cstdlib>
#include <iostream>
#include <math.h>
#include "mersenne.cpp"             // code for random number generator
#include "userintf.cpp"             // define system specific user interface
#include <fstream>
#include <iomanip>
#include "equation_SFH.c"
#include "SNII_distance.c"
#include "SNIa_distance.c"
#include <string>
#include <sstream>
#include <stdlib.h>
#include "SFH/SFH.cpp"

using namespace std;

//global constants
/* ******************************************** */  
  
  //subsection slice distance in terms of radial distance
  #define subsection_slice_radial_distance 30; //in parsecs
  
  //subsection slice distance in terms of height
  #define subsection_slice_h_z_distance 30; //in parsecs
  
  //radius of the milky way disk
  #define mw_radius 15330; //in parsecs

  //height of the milky way disk (total thickness)
  //5000 captures 99% of all stars in a 10000 pc vertical height
  //#define mw_height 5000; //in parsecs
  #define mw_height 5000; //in parsecs- for testing
  
  //angle
  #define mw_angle 1; //in degrees
  
  //PI
  #define PI 3.14159265;
  
  //*************
  
  #define number_density_type 1; //1-C&O
  //stellar density*************************
  //in parsecs- Carroll and Ostlie
  #define h_R 2250;
  #define Z_thin 350;
  #define Z_thick 1000;
  
  // #define n0 5.502; //for testing
  // #define n0 5.502; //this is the appropriate value corresponding to Salpeter IMF



  //stellar mass-Salpeter
  #define m_max 100; //maximum stellar mass allowed in solar masses
  #define m_min 0.08;//minimum stellar mass allowed in solar masses
  #define alpha 2.35; //exponent in the salpeter IMF
  
 
//*******************************  

  //SN sterilize distance
  #define typeII 8; //in parsecs
  #define typeIa 20;
  //mass that causes a typeII SN in solar masses
  #define typeII_mass 8; 
  
  //lower and upper limits of stellar masses that could cause a Type Ia SN
  #define type_Ia_mass_lower 0.9285427; //this corresponds to 13.24 Byr: if a star lives longer than this period it will not become a SNIa in the simulation
                                       //eg: an M class star that will live to be 1000+ billion years is not a threat
  #define type_Ia_mass_upper 8;

  //time for complex life to evolve************
  #define max_time 13.24; //the temporal length (oldest known star in the MWG is 13.2 years old
                          //the .04 extra is due to rounding errors in the SFH
                          
  
  //sterilization length
  #define max_len_ster_hist 10; //10 elements is the maximum number that can be stored
  
  //sterilization distance multiplier
  #define ster_dist_multi 1; //keep this at 1- corresponding to an 8pc normalization to SNII

struct subsection_parameters
{
       //which of the xy plane subsections it belongs to
       int subsection_xy; 
       //which of the z plane the subsection belongs to
       int subsection_z;
       //to the centre of the subsection
       float mp_z_coord_height_above_midplane;
       //to the centre of the subsection
       float midpoint_xy_plane_dist;
       //volume of the subsection
       float volume;

       float radial_distance_min;
       float radial_distance_max;
       float angle_min; //this is always 0- included for descriptive purposes
       float angle_max;
       float z_min;
       float z_max;     
       
       long int num_stars;
       
       long int star_index_min;
       long int star_index_max;
};



struct star_cell
{
       int subsection_z;
       int cell;   
       float x_coord;
       float y_coord;   
       float z_coord;      
       int sn;     
       int sterilized;
       int sterilized_count; //number of times a star is sterilized              
       float mass;
       float birth_date;
       float ms_lifetime;
       float sterilization_date;
       float death_date; 
       float sterilization_distance; //only applicable to stars that are SN-holds the sterilization distance for a SNIa or SNII
};

struct column_properties
{
       int cell_min;
       int cell_max;
       long int num_stars; //for dynamic memory allocation
       
};

//stats structure
struct statistics
{
            long long total_sterilizations;   
            long long sterilizations_typeII;   
            long long sterilizations_typeIa;   
            float distance;    
            long long total_stars;
            long long habitable; //sterilized but habitable for 4gyr
            long long not_sterilized_and_habitable;

          
};



/* ******************************************** */  
//global variables

  //pointer to a dynamic array of structures containing the parameters for each subsection
  subsection_parameters **subsect_ptr;
  int count_z_glbl;
  int count_xy_glbl;
  unsigned long long num_stars_glbl=0;
  unsigned long long num_sn_glbl=0;
  unsigned long long num_sn_Ia_glbl=0;
  
  unsigned long long num_SNII_specific=0;
  unsigned long long num_SNIa_specific=0;
  long column_stars=0;
  int stats_counter=0;
  long double total_stellar_mass_glbl=0;
  int ptr_1_range;
  int ptr_2_range;
  int ptr_3_range;
  int delete_ptr;
  int star_cell_ptr_1, star_cell_ptr_2, star_cell_ptr_3, delete_ptr_range;
  
star_cell *star; 
star_cell *star1;
star_cell *star2;
star_cell wrap_around; //used for the wrap around- periodic boundary effect, only need 1 structure
  
  //outlines the properties of each cell
  struct subsection_parameters **subsect_param_array = NULL;
  
  //the number of stars per column
  struct column_properties **column_array=NULL;
  
  unsigned long int total_sterilizations_II=0;
  unsigned long int total_sterilizations_Ia=0;
  
  
 //stats
 //c allocation
 struct statistics **stats = NULL;
  
//c++ allocation  
statistics *stats_midplane;
statistics *stats_radial_mp;// stats for the number of habitable planets as a function of radial distance and height above the midplane
                               //this needs its own stats because radial and midplane stats are implemented together, not seperately

//stats for the galactic evolution
long double num_stars_1byr_increments[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};//holds the number of stellar mass for particular time periods cumulative: 0Bya (everything that has formed), 1Bya -1 billion years ago to the formation of the Galaxy etc

//function prototypes
/* ******************************************** */  
int determine_habitable_birthdate(int index);
void create_subsections();
float random_star_generator_R (float deg_min, float deg_max, float r_min, float r_max);
float random_star_generator_Z(float z_min, float z_max);
float degree_generator(float deg_min, float deg_max);
int SN_star(float mass); 
int wrap_around_function(float x, float y, float z, float sterilization_distance);
long int get_star_range_lower(int sn_star);
long int get_star_range_upper(int sn_star);
long int get_star_range_lower_cell_right(int cell);
long int get_star_range_upper_cell_right(int cell);
long int get_star_range_lower_cell_left(int cell);
long int get_star_range_upper_cell_left(int cell);
long int whole_num_star(float number);
void print_stats();//outputs radial stats to the stats structure
void output_stats(); //outputs stats to the stats file
float stellar_mass();
float Kroupa_IMF();
float Salpeter_IMF();
void create_columns ();
void copy_into_star1(long int length);
void copy_into_star2(long int length);
void copy_star2_into_star1(long int length);
void copy_temp_into_star2(long int length);
void process_SN();
void process_SN_end_subsection();
long int running_total(int column);//returns the number of stars before the column
float generate_ms_lifetime(float mass);//determines the lifetime of a star
int type_Ia_candidate();//determines if a star that will become a white dwarf will also be a type Ia SN
//int evaluate_sterilization_history_planet(int index, float hab_time); 
//float evaluate_when_planet_habitable(int index);
int CmpFunc(const void* _a, const void* _b);//for quicksort

//random number generation must be global- seeding the random number generator multiple times produces the same 
//random numbers
//int32 seed = (int32)time(0);        // random seed
int32 seed = (int32)666;        // random seed
CRandomMersenne rg(seed);           // make instance of random number generator
/* ******************************************* */

// Global variable
float n0;

int main(int argc, char *argv[])
{
    
   
    int state;
    int num_in_process=0;
    int num_iterated=0;
    unsigned long int count=1;
    float R; //radial distance
    float degree; //dandom degree
    float pi_fn=PI;
    int half_cells_for_stats=0;
    
    float type_Ia_mass_upper_fn =type_Ia_mass_upper;
    float type_Ia_mass_lower_fn =type_Ia_mass_lower;
    float max_time_fn=max_time; 
    float radial_dist_cell_size=subsection_slice_radial_distance;
    float sterilization_multiplier=ster_dist_multi;

    if (argc != 2) {
      fprintf(stderr,"Usage: %s <n0>\n",argv[0]);
      fprintf(stderr,"    The larger n0 the longer the runtime (default n0 = 5.502)\n");
      exit(1);
    }
    if ((sscanf(argv[1], "%f", &n0) != 1) || (n0 < 0))  {
      fprintf(stderr,"    Invalid command-line argument (should be a positive float)\n");
      exit(1);
    }

    create_subsections();
    create_columns();
    

    
    cout <<"\n count xy: "<<count_xy_glbl;
    cout <<"\n count z: "<<count_z_glbl;
    cout << "\n number stars: "<<num_stars_glbl;
    cout <<"\n";
    
    for (int l=1; l<=count_xy_glbl; l++)
    {

            
            star = new star_cell [column_array[l]->num_stars+1];
            //loop that iterates through each cell in the column
            for (long int j=column_array[l]->cell_min; j<=column_array[l]->cell_max; j++)
            {
                
                    for (int k=1; k<=subsect_ptr[j]->num_stars; k++)
                    {
                
                    star[count].subsection_z=subsect_ptr[j]->subsection_z;
                    
                    star[count].cell=j;
                    //for radial distance in cylindrical coordinates
                    R=random_star_generator_R(subsect_ptr[j]->angle_min,subsect_ptr[j]->angle_max, 
                    subsect_ptr[j]->radial_distance_min, subsect_ptr[j]->radial_distance_max);
                    
                    degree=degree_generator(subsect_ptr[j]->angle_min, subsect_ptr[j]->angle_max);
                    //set x,y,z coordinates
                    star[count].x_coord=R*(cos(degree*(pi_fn/180)));
                    star[count].y_coord=R*(sin(degree*(pi_fn/180)));
                    star[count].z_coord=random_star_generator_Z(subsect_ptr[j]->z_min,subsect_ptr[j]->z_max);
                    
                    //initialize sterilization value to 0
                    star[count].sterilized=0;
                    //initialize sn value to 0
                    star[count].sn=0;
                    star[count].mass=stellar_mass();
                    total_stellar_mass_glbl+=star[count].mass;

                    star[count].birth_date=generate_birth_date((R/1000),rg.Random());            
                    star[count].ms_lifetime=generate_ms_lifetime(star[count].mass);
                    star[count].death_date=star[count].birth_date+star[count].ms_lifetime;
                                    
                    //determines if a star is a typeII SN, returns 1 if it is a SN
                    if (
                    (SN_star(star[count].mass))==1
                    )
                    {
                    star[count].sn=1; 
                    //sets the sterilization distance of the SNII
                    star[count].sterilization_distance=SNII_distance(rg.Random());
                    
                    //multiplies the serilization distance by the multiplier for sensitivity analysis
                    star[count].sterilization_distance*=sterilization_multiplier;
                    
                    if (star[count].death_date<=13.24)
                          {  
                          num_sn_glbl++; 
                          
                            //calculate some global statistics between 2.5-15kpc
                            //to remove those statistics between 0-2.5kpc
                               if (l>(int((2500/radial_dist_cell_size))) )
                               {
                                   num_SNII_specific++; 
                               }  
                          }      
                      
                  
                        
                    }
                    
                    //determines if a star will cause a typeIa SN
                    //it will return 2 if it is, and 0 if it isn't
                    if (
                    (star[count].mass >=type_Ia_mass_lower_fn)&&(star[count].mass <type_Ia_mass_upper_fn) && (star[count].death_date<=max_time_fn)
                    )
                    {
                    star[count].sn=type_Ia_candidate();       
                    
                      //if the star is a typeIa SN- set its sterilization distance
                      if (star[count].sn==2)
                      {
                         //sets the sterilization distance of the SNIa
                         star[count].sterilization_distance=SNIa_distance(rg.Random());   
                          //multiplies the serilization distance by the multiplier for sensitivity analysis
                         star[count].sterilization_distance*=sterilization_multiplier;
                         
                          num_sn_Ia_glbl++;    
                          
                          //calculate some global statistics between 2.5-15kpc
                          //to remove those statistics between 0-2.5kpc
                          if (l>(int((2500/radial_dist_cell_size))) )
                          {
                          num_SNIa_specific++;   
                                  
                          }                          

                                       
                      }                                                
                    }
              
                    
                    count++;     
                                          
                       }
            
            
           
        } //end of loop
        
         num_in_process++;
        
        
        
        //special case where the first column needs to be populated
        if (num_in_process==1)
        {
           ptr_1_range=count-1; 
           copy_into_star1(ptr_1_range);           
          
          //initialize stats array
          stats = (struct statistics **)realloc(stats, (count_xy_glbl+1) * sizeof(struct statistics *));
          

            stats_midplane = new statistics [count_z_glbl+1];  
            stats_radial_mp = new statistics [(count_z_glbl*count_xy_glbl)+1];
                                                     
        }
        //need num_in_process2 to populate star2 with star
        else if (num_in_process==2)
        
        {
            ptr_2_range=count-1;
            copy_into_star2(ptr_2_range);
            cout <<"\nnum in process: "<<num_in_process<<"/"<<count_xy_glbl;;    
            process_SN();     
            print_stats();  
          //  print_stats_midplane(); 
        }
        
        //if it is the last column being processed
        //do: 1->1, 1->2, 2->1, 2->2
        else if (num_in_process==count_xy_glbl)
        {
             
               //like normal
                ptr_1_range=ptr_2_range;
                ptr_2_range=count-1; 
                
                copy_star2_into_star1(ptr_1_range); 
                copy_temp_into_star2(ptr_2_range); 
                cout <<"\nnum in process: "<<num_in_process<<"/"<<count_xy_glbl;
                process_SN();    
                print_stats(); 
              //  print_stats_midplane(); 
                //now run the process SN on itself 
                
               
                copy_star2_into_star1(ptr_2_range);
                ptr_1_range=ptr_2_range;
                cout <<"\nnum in process: "<<num_in_process<<"/"<<count_xy_glbl;;
                process_SN_end_subsection();
                print_stats(); 
             //   print_stats_midplane(); 

                //output will show num in process 2 times, cause the last one does get processed twice
             
        }
        
        //this is the normal condition
            else
            {
                //shift the dynamic array star2 to star1 and the new column into star2
                ptr_1_range=ptr_2_range;
                ptr_2_range=count-1; 
                delete [] star1;
                copy_star2_into_star1(ptr_1_range);
                delete [] star2;
                copy_temp_into_star2(ptr_2_range);
                delete [] star;
                //copy_into_star2(ptr_2_range); 
                cout <<"\nnum in process: "<<num_in_process<<"/"<<count_xy_glbl;
                process_SN();    
                print_stats(); 

      
            }

         count=1;
  
    


}//end of loop


    
   
cout << "\n number stars: "<<num_stars_glbl;
cout << "\n number of SNII: "<<num_sn_glbl;
cout << "\n number of SNIa: "<<num_sn_Ia_glbl;
cout << "\n total stellar mass: "<<total_stellar_mass_glbl;
cout << "\n Sterilization distance multiplier: "<<ster_dist_multi;
cout << "\n Number density type: "<<number_density_type;
cout << "\n";


    output_stats();
  
    //remember to delete the dynamic arrays columns
     delete subsect_ptr;
     delete star,star1,star2;
     delete column_array;
     delete stats;
     
    return 0;
}

void process_SN_end_subsection()
{
     float sn_x, sn_y, sn_z;
     
     float distance;
     float distance2;
     float sn2=typeII;
     float snIa=typeIa;
     int sn_type=0;
     int wrap_flag;
     float sterilization_distance;
     float birthdate_sn=0;
     float birthdate_candidate=0;
     float deathdate_sn=0;
     float deathdate_candidate=0;
     float max_time_fn=max_time;
      int normal_ster_flag=0;
      int wrap_ster_flag=0;
     
     //always looking for the SN as they occur in star1
     //eg 1->1, 1->2, 2->1
     
         //********************************
         //1->1 on itself and 1->2
         for (long int i=1; i<=ptr_1_range; i++)
         {
             
             //********************************
         //1->1 on itself
         //*************************************
                     if(star1[i].sn==1 || star1[i].sn==2)
                     {
                              
                             //set the birthdate/deathdate of the SN star
                             birthdate_sn=star1[i].birth_date;
                             deathdate_sn=star1[i].death_date;
                               //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             if (star1[i].sn==1)
                             {           
                             sterilization_distance=star1[i].sterilization_distance; 
                             sn_type=1;
                             }          
                             //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             else if (star1[i].sn==2)
                             {
                             sterilization_distance=star1[i].sterilization_distance;   
                             sn_type=2;  
                             }
                                        
                               
                              sn_x=star1[i].x_coord;
                              sn_y=star1[i].y_coord;
                              sn_z=star1[i].z_coord;
                              wrap_flag=0;
                              if(wrap_around_function(sn_x, sn_y, sn_z,sterilization_distance)==1)
                              {
                              wrap_flag=1;                              
                              }
                                 
                              for (long int j=get_star_range_lower(i); j<=get_star_range_upper(i); j++)
                              {
                                  
                                  normal_ster_flag=0;
                                  wrap_ster_flag=0;
                                  birthdate_candidate=star1[j].birth_date;
                                  deathdate_candidate=star1[j].death_date;
                                  distance=(sqrt
                                               (
                                               ((sn_x-star1[j].x_coord)*(sn_x-star1[j].x_coord))+
                                               ((sn_y-star1[j].y_coord)*(sn_y-star1[j].y_coord))+
                                               ((sn_z-star1[j].z_coord)*(sn_z-star1[j].z_coord))
                                               ));
                                  
                                  //wrap around distance
                                  if (wrap_flag==1)
                                  {
                                  distance2=(sqrt
                                               (
                                               ((wrap_around.x_coord-star1[j].x_coord)*(wrap_around.x_coord-star1[j].x_coord))+
                                               ((wrap_around.y_coord-star1[j].y_coord)*(wrap_around.y_coord-star1[j].y_coord))+
                                               ((wrap_around.z_coord-star1[j].z_coord)*(wrap_around.z_coord-star1[j].z_coord))
                                               ));                               
                                                                 
                                  }
                                  
                                  //normal sterilization
                                  if(
                                  (distance<=sterilization_distance) && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn)
                                  && (deathdate_sn <=max_time_fn)
                                  )
                                  {
                                      if (sn_type==1)
                                     {
                                     star1[j].sterilized=1;  
                                     }
                                     else if (sn_type==2)
                                     {
                                     star1[j].sterilized=2;  
                                     }
 
                                     star1[j].sterilization_date=deathdate_sn; 
                                     normal_ster_flag=1;       
                                      
                                     
                                
                                     
                                      
                                  }
                                    
                                    //wrap around sterilization
                                    if(
                                    (wrap_flag==1) && (distance2<=sterilization_distance)
                                     && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn) 
                                     && (deathdate_sn <=max_time_fn)
                                     )
                                  {
                                     if (sn_type==1)
                                     {
                                     star1[j].sterilized=1;  
                                     }
                                     else if (sn_type==2)
                                     {
                                     star1[j].sterilized=2;  
                                     }

                                     star1[j].sterilization_date=deathdate_sn;
                                     
                                     wrap_ster_flag=1;
                                  
                                     
                                      
                                  } // end of wrap around sterilization
                                  
                                  if(normal_ster_flag==1 || wrap_ster_flag==1)
                                    {
                                          
                                          if (sn_type==1)
                                          {
   
                                           star1[j].sterilized_count++;
                                           total_sterilizations_II++;
                                           }           
                                           else if (sn_type==2)
                                           {    
                                           total_sterilizations_Ia++;

                                           } 
                                    }  

                                  
                                  
                              } //end of loop
                     }
         }
}

void process_SN()
{ 
     float sn_x, sn_y, sn_z;
     int wrap_flag=0;
     float distance;
     float distance2;//distance between SN wrap around star, and the star being checked
     float sn2=typeII;
     float snIa=typeIa;
     int sn_type=0;
     float sterilization_distance;
     float birthdate_sn=0;
     float birthdate_candidate=0;
     float deathdate_sn=0;
     float deathdate_candidate=0;
     float max_time_fn=max_time;
     int normal_ster_flag=0;
     int wrap_ster_flag=0;
     //always looking for the SN as they occur in star1
     //eg 1->1, 1->2, 2->1
     
         //********************************
         //1->1 on itself and 1->2
         for (long int i=1; i<=ptr_1_range; i++)
         {
            
             //********************************
         //1->1 on itself
         //*************************************
                     if(star1[i].sn==1 ||star1[i].sn==2)
                     {
                            
                             //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             if (star1[i].sn==1)
                             {           
                             sterilization_distance=star1[i].sterilization_distance; 
                             sn_type=1;
                             }          
                             //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             else if (star1[i].sn==2)
                             {
                             sterilization_distance=star1[i].sterilization_distance;   
                             sn_type=2;  
                             }
                             
                             //set the birthdate/deathdate of the SN star
                             birthdate_sn=star1[i].birth_date;
                             deathdate_sn=star1[i].death_date;
                            
                              sn_x=star1[i].x_coord;
                              sn_y=star1[i].y_coord;
                              sn_z=star1[i].z_coord;
                              
                              
                              wrap_flag=0;

                              if(wrap_around_function(sn_x, sn_y, sn_z,sterilization_distance)==1)
                              {
                              wrap_flag=1;                              
                              }
                              
                                  //  sn_debug << "\nstar range lower self "<<get_star_range_lower(i);
                                   // sn_debug << "\nstar range upper self: "<<get_star_range_upper(i);
                              for (long int j=get_star_range_lower(i); j<=get_star_range_upper(i); j++)
                              {
                                  normal_ster_flag=0;
                                  wrap_ster_flag=0;
                                  birthdate_candidate=star1[j].birth_date;
                                  deathdate_candidate=star1[j].death_date;
                                  
                                  //normal condition-calculates distances between stars within the cell
                                  distance=(sqrt
                                               (
                                               ((sn_x-star1[j].x_coord)*(sn_x-star1[j].x_coord))+
                                               ((sn_y-star1[j].y_coord)*(sn_y-star1[j].y_coord))+
                                               ((sn_z-star1[j].z_coord)*(sn_z-star1[j].z_coord))
                                               ));
                                  //SN star is within 10pc from the x-axis or the (1) degree incline
                                  //wrap around
                                  if (wrap_flag==1)
                                  {
                                  distance2=(sqrt
                                               (
                                               ((wrap_around.x_coord-star1[j].x_coord)*(wrap_around.x_coord-star1[j].x_coord))+
                                               ((wrap_around.y_coord-star1[j].y_coord)*(wrap_around.y_coord-star1[j].y_coord))+
                                               ((wrap_around.z_coord-star1[j].z_coord)*(wrap_around.z_coord-star1[j].z_coord))
                                               ));                               
                                                                 
                                  }
                                  
                                  
                                  //normal sterilization
                                  if(
                                  (distance<=sterilization_distance) 
                                  && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn)  
                                  && (deathdate_sn <=max_time_fn)
                                  
                                  )
                                  {
                                     
                                     if (sn_type==1)
                                     {
                                     star1[j].sterilized=1;  
                                     
                                     }
                                     else if (sn_type==2)
                                     {
                                     star1[j].sterilized=2;  
                                    
                                     }
                                      
                                      
                                     star1[j].sterilization_date=deathdate_sn;
                                     normal_ster_flag=1;
                                     
                                     
                                    
                                  } //end of normal sterilization
                                  
                                   
                                  
                                  //wrap around sterilization
                                  if (
                                  (wrap_flag==1) && (distance2<=sterilization_distance)&& 
                                  (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn) 
                                  && (deathdate_sn <=max_time_fn)
                                  
                                  )
                                  {
                                        if (sn_type==1)
                                         {
                                         star1[j].sterilized=1;
                                          
                                         }
                                         else if (sn_type==2)
                                         {
                                         star1[j].sterilized=2;  
                                        
                                         }

                                         star1[j].sterilization_date=deathdate_sn;
                                         wrap_ster_flag=1;
                                   
                                                   
                                  } //end of wrap sterilization
                                  
                              
                                  
                                   if(normal_ster_flag==1 || wrap_ster_flag==1)
                                    {
                                          
                                          if (sn_type==1)
                                          {  
                                           star1[j].sterilized_count++;
                                           total_sterilizations_II++;
                                           }           
                                           else if (sn_type==2)
                                           {
                                            star1[j].sterilized_count++; 
                                            total_sterilizations_Ia++;
                                           } 
                                    }  
                                      
                                       
                              
                                  
                              } //end of loop for sterilizations 1->1
                              
                              
         //********************************
         //1->2 
         //*************************************                      
                          //   sn_debug << "\n1->2 "<<endl; 
                            // sn_debug << "star: "<<i<<" went SN"<<"  Cell: "<<star1[i].cell<<"  Subsection: "<<star1[i].subsection_z<<"  Column: "<<subsect_param_array[star1[i].cell]->subsection_xy;   
                                //        sn_debug << "\nstar range lower right: "<<get_star_range_lower_cell_right(star1[i].cell);
                               //     sn_debug << "\nstar range upper right: "<<get_star_range_upper_cell_right(star1[i].cell);
                              for (long int k=get_star_range_lower_cell_right(star1[i].cell); k<=get_star_range_upper_cell_right(star1[i].cell); k++)
                              {
                                  normal_ster_flag=0;
                                  wrap_ster_flag=0;
                        
                                  birthdate_candidate=star2[k].birth_date;
                                  deathdate_candidate=star2[k].death_date;
                                  //normal condition
                                  distance=(sqrt
                                               (
                                               ((sn_x-star2[k].x_coord)*(sn_x-star2[k].x_coord))+
                                               ((sn_y-star2[k].y_coord)*(sn_y-star2[k].y_coord))+
                                               ((sn_z-star2[k].z_coord)*(sn_z-star2[k].z_coord))
                                               ));
                                               
                                    
                                  //if the star wraps around             
                                  if (wrap_flag==1)
                                  {
                                  distance2=(sqrt
                                               (
                                               ((wrap_around.x_coord-star2[k].x_coord)*(wrap_around.x_coord-star2[k].x_coord))+
                                               ((wrap_around.y_coord-star2[k].y_coord)*(wrap_around.y_coord-star2[k].y_coord))+
                                               ((wrap_around.z_coord-star2[k].z_coord)*(wrap_around.z_coord-star2[k].z_coord))
                                               ));                               
                                                                 
                                  }
                                  
                                  //normal sterilization
                                  if(
                                  (distance<=sterilization_distance) && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn)
                                  && (deathdate_sn <=max_time_fn)
                                  )
                                  {
                                  
                                        if (sn_type==1)
                                         {
                                         star2[k].sterilized=1;  
                                         }
                                         else if (sn_type==2)
                                         {
                                         star2[k].sterilized=2;  
                                         }
                                         star2[k].sterilization_date=deathdate_sn;
                                         normal_ster_flag=1;
                                           //for debugging************
                                           
                                
                                   //  sn_debug << "\nsterilization at star2: "<<k;
                                   //  sn_debug <<"\nCell: "<<star2[k].cell<<"  Subsection: "<<star2[k].subsection_z<<"  Column: "<<subsect_param_array[star2[k].cell]->subsection_xy;
                                   //  sn_debug <<"\n1->2 sterilized by star1: "<<i<<" cell: "<<star1[i].cell;
                                    // sn_debug << "\nsterilization(normal) at star2 time: "<<star2[k].sterilization_date<<endl<<endl;
                                     /*
                                     sn_debug << "\nsterilization map: ";
                                     for (int h=1; h<=133; h++)
                                     {
                                     sn_debug <<star2[k].ster_map[h];
                                     }
                                     sn_debug<<endl;
                                     */
                                    //**************************************** 
                                    
                                       
                                           
                                  }
                                  
                                 
                                 
                                  //wrap-around sterilization
                                  if(
                                  (wrap_flag==1) && (distance2<=sterilization_distance)
                                  && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn) 
                                  && (deathdate_sn <=max_time_fn)
                                  )
                                  {
                                  
                                         if (sn_type==1)
                                         {
                                         star2[k].sterilized=1;  
                                        
                                         }
                                         else if (sn_type==2)
                                         {
                                         star2[k].sterilized=2;  
                                         
                                         }
                                         star2[k].sterilization_date=deathdate_sn;
                                         wrap_ster_flag=1;

                                      
                                             
                                      
                                  } //end of wrap sterilization
                                  
                                  
                                 
                                              if(normal_ster_flag==1 || wrap_ster_flag==1)
                                              {
                                                  
                                                  if (sn_type==1)
                                                  {
  
                                                   star2[k].sterilized_count++;
                                                 
                                                   total_sterilizations_II++;
                                                   }           
                                                   else if (sn_type==2)
                                                   {

                                                    star2[k].sterilized_count++;  
                                                    total_sterilizations_Ia++;

                                                   } 
                                              }  
                                           
                              }  //end of for loop for 1->2 sterilizations
                            
                                        
                     } //end of if sterilization in star1
                     
          
          
       
         
         }   //end of 1->1 and 1->2- the loop that scans through all of star1 
         
           //********************************
         //2->1 
         //*************************************   
         
         for (long int l=1; l<=ptr_2_range; l++)
         {
             
                  if(star2[l].sn==1 || star2[l].sn==2)
                  {
                              
                                
                             //set the birthdate/deathdate of the SN star
                             birthdate_sn=star2[l].birth_date;
                             deathdate_sn=star2[l].death_date;
                             //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             if (star2[l].sn==1)
                             {           
                             sterilization_distance=star2[l].sterilization_distance; 
                             sn_type=1;
                             }          
                             //sterilization distance for typeII
                             //to be replaced with the distributions of distances
                             else if (star2[l].sn==2)
                             {
                             sterilization_distance=star2[l].sterilization_distance;   
                             sn_type=2;  
                             }       
                                   
                       sn_x=star2[l].x_coord;
                       sn_y=star2[l].y_coord;
                       sn_z=star2[l].z_coord;  
                       
                              wrap_flag=0;
                              if(wrap_around_function(sn_x, sn_y, sn_z, sterilization_distance)==1)
                              {
                              wrap_flag=1;                              
                              }
                       
                
                           for (long int m=get_star_range_lower_cell_left(star2[l].cell); m<=get_star_range_upper_cell_left(star2[l].cell); m++)
                              {
                                  normal_ster_flag=0;
                                  wrap_ster_flag=0;    
                                  birthdate_candidate=star1[m].birth_date;
                                  deathdate_candidate=star1[m].death_date;
                                  //normal sterilization distance
                                  distance=(sqrt
                                               (
                                               ((sn_x-star1[m].x_coord)*(sn_x-star1[m].x_coord))+
                                               ((sn_y-star1[m].y_coord)*(sn_y-star1[m].y_coord))+
                                               ((sn_z-star1[m].z_coord)*(sn_z-star1[m].z_coord))
                                               ));
                                         
                                   //wrap around sterilization distance            
                                   if (wrap_flag==1)
                                  {
                                  distance2=(sqrt
                                               (
                                               ((wrap_around.x_coord-star1[m].x_coord)*(wrap_around.x_coord-star1[m].x_coord))+
                                               ((wrap_around.y_coord-star1[m].y_coord)*(wrap_around.y_coord-star1[m].y_coord))+
                                               ((wrap_around.z_coord-star1[m].z_coord)*(wrap_around.z_coord-star1[m].z_coord))
                                               ));                               
                                                                 
                                  }
                                  
                                  //normal sterilization
                                  if(
                                   (distance<=sterilization_distance) && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn)
                                   && (deathdate_sn <=max_time_fn)
                                   )
                                  {
                                        if (sn_type==1)
                                         {
                                         star1[m].sterilized=1;  
                                         }
                                         else if (sn_type==2)
                                         {
                                         star1[m].sterilized=2;   
                                         }
                                      star1[m].sterilization_date=deathdate_sn;
                                      normal_ster_flag=1;
                                      
                                     
                                     
                                    
                                     
                                  }
                                  
                                
                                 
                                  //wrap around sterilization
                                   if(
                                   (wrap_flag==1) && (distance2<=sterilization_distance)
                                   && (birthdate_candidate<deathdate_sn) && (deathdate_candidate > deathdate_sn)
                                    && (deathdate_sn <=max_time_fn)
                                   )
                                   {
                                       if (sn_type==1)
                                         {
                                        star1[m].sterilized=1;   
                                         }
                                         else if (sn_type==2)
                                         {
                                         star1[m].sterilized=2;   
                                         }
                                      star1[m].sterilization_date=deathdate_sn;
                                      wrap_ster_flag=1;
                                      
                        
                                     
                                    //****************************************  
                                    
                                    
                                                    
                                                   
                                     } //end of wrap around sterilization
                                  
                                  
                                 
                                 
                                                if(normal_ster_flag==1 || wrap_ster_flag==1)
                                                  {
                                                  
                                                  if (sn_type==1)
                                                  {
                                                   star1[m].sterilized_count++;
                                                  // update_sterilization_history(m,0,deathdate_sn);
                                                   //update_SN_sterilization_history(m,0,1);
                                                   total_sterilizations_II++;
                                                   }           
                                                   else if (sn_type==2)
                                                   {   
                                                     star1[m].sterilized_count++;    
                                                     //update_sterilization_history(m,0,deathdate_sn);
                                                     //update_SN_sterilization_history(m,0,2);
                                                     total_sterilizations_Ia++;
                                                   } 
                                                   } 
                                                   
                                                      //*********
                                
                              } //end of loop that scans through looking to see if a star gets sterilized in star1 2->1
                       
                                        
                  } //end of if statement that checks to see if a star was sterilized
                                     
                                     
         } //end of loop that scans through all of the stars to see if they are sterilized
         
}

int type_Ia_candidate()
{
    
    if (rg.Random()<=0.01)
    {
     
     return 2;
       
    }
    else
    {
        return 0;
    }
    
    
}

//determines if a star is close to a border, if it is return 1
//change the coordinates in star wrap_around and perform the appropriate actions in process_SN()
int wrap_around_function(float x, float y, float z, float sterilization_distance)
{
    float R;
    float angle;
    float PI_fn=PI;
    float arc_length_original;
    float arc_length_total;
    float arc_length_ratio;
    float mw_angle_fn=mw_angle;
    float outside_angle;

    
    R=sqrt(
    (x*x)+(y*y)
    );
    
    angle=(acos(x/R))*(180/PI_fn);
    
    arc_length_original=(angle*(2*PI_fn/360))*(R);
    arc_length_total=(mw_angle_fn*(2*PI_fn/360))*(R); 
    arc_length_ratio=arc_length_original/arc_length_total;
    
    
    //determine if the SN star passed to this function is within the distance passed to this function to one of the borders
    //of either boundary
    //if it isn't return 0
    if (
    !((arc_length_total-arc_length_original)<sterilization_distance) && !(arc_length_original<sterilization_distance)
    )
    {
           //  cout <<"\nhere1!";                        
       return 0;
    }
                                    
    //the SN star is closest to the x-axis
    //wraps around to the outside 1 degree incline line
    //shouldnt be <=, however, there's virtually no chance it will ever be exactly 0.5
    if (arc_length_ratio<=0.5)
    {
        outside_angle=mw_angle_fn+(mw_angle_fn*angle);    
        
        wrap_around.x_coord=R*(cos(outside_angle*(PI_fn/180)));                
        wrap_around.y_coord=R*(sin(outside_angle*(PI_fn/180))); 
        wrap_around.z_coord=z;
         return 1;                     
    }
    //the SN star is closest to the 1 degree incline line
    //wraps around to the x-axis
    else //(arc_length_ratio>0.5)
    {
        outside_angle=mw_angle_fn-angle;  //below the x-axis, from the origin
      
        wrap_around.x_coord=R*(cos(outside_angle*(PI_fn/180)));
        //make the y coordinate negative so it's below the x-axis  
        wrap_around.y_coord=-1*(R*(sin(outside_angle*(PI_fn/180))));    
        wrap_around.z_coord=z;        
         return 1;   
    }    
          
}

float generate_ms_lifetime(float mass)
{
int sun_ms_lifetime=11;
int m_sun=1;
float ms_lifetime;

ms_lifetime= (pow((m_sun/mass),2.5))*sun_ms_lifetime;

return(ms_lifetime);
      
}
//generates the time a star formed
float generate_SFH()
{
      float current_time=max_time; //13.24
double x;
x=rg.Random();
return(current_time-eqn7909(x));      
}

//returns the total number of stars before the column
long int running_total(int column)
{
     long int total=0;
     
     for (int i=1; i<column; i++)
     {
            total=total+column_array[i]->num_stars; 
     }
     
     return (total);

}

void copy_into_star1(long int length)
{
     int max_len_ster_hist_fn=max_len_ster_hist;
    star1 = new star_cell [length+1];
     for (long int i=1; i<=length; i++)
     {
         
         star1[i].cell=star[i].cell;    
         star1[i].mass=star[i].mass; 
         star1[i].sn=star[i].sn; 
         star1[i].sterilized=star[i].sterilized; 
         star1[i].subsection_z=star[i].subsection_z; 
         star1[i].x_coord=star[i].x_coord; 
         star1[i].y_coord=star[i].y_coord; 
         star1[i].z_coord=star[i].z_coord;       
         star1[i].birth_date=star[i].birth_date;
         star1[i].death_date=star[i].death_date;
         star1[i].ms_lifetime=star[i].ms_lifetime;
         star1[i].sterilization_date=star[i].sterilization_date;
         star1[i].sterilized_count=star[i].sterilized_count;
         star1[i].sterilization_distance=star[i].sterilization_distance;
         
     }
     
     
     
     
}


void copy_into_star2(long int length)
{
int max_len_ster_hist_fn=max_len_ster_hist;
    star2 = new star_cell [length+1];
     for (long int i=1; i<=length; i++)
     {

         
         star2[i].cell=star[i].cell;    
         star2[i].mass=star[i].mass; 
         star2[i].sn=star[i].sn; 
         star2[i].sterilized=star[i].sterilized; 
         star2[i].subsection_z=star[i].subsection_z; 
         star2[i].x_coord=star[i].x_coord; 
         star2[i].y_coord=star[i].y_coord; 
         star2[i].z_coord=star[i].z_coord;
         
         star2[i].birth_date=star[i].birth_date;
         star2[i].death_date=star[i].death_date;
         star2[i].ms_lifetime=star[i].ms_lifetime;
         star2[i].sterilization_date=star[i].sterilization_date;
         star2[i].sterilized_count=star[i].sterilized_count;
         star2[i].sterilization_distance=star[i].sterilization_distance;
     }
     
    
}


void copy_star2_into_star1(long int length)
{
int max_len_ster_hist_fn=max_len_ster_hist;



star1 = new star_cell [length+1];
     for (long int i=1; i<=length; i++)
     {
   
         
         star1[i].cell=star2[i].cell;    
         star1[i].mass=star2[i].mass; 
         star1[i].sn=star2[i].sn; 
         star1[i].sterilized=star2[i].sterilized; 
         star1[i].subsection_z=star2[i].subsection_z; 
         star1[i].x_coord=star2[i].x_coord; 
         star1[i].y_coord=star2[i].y_coord; 
         star1[i].z_coord=star2[i].z_coord;
         
         star1[i].birth_date=star2[i].birth_date;
         star1[i].death_date=star2[i].death_date;
         star1[i].ms_lifetime=star2[i].ms_lifetime;
         star1[i].sterilization_date=star2[i].sterilization_date;
         star1[i].sterilized_count=star2[i].sterilized_count;
         star1[i].sterilization_distance=star2[i].sterilization_distance;
     
     } 
      

}

void copy_temp_into_star2(long int length)
{
int max_len_ster_hist_fn=max_len_ster_hist;
//*star2 = NULL;

//star2 = (struct star_cell **)realloc(star2, (length+1) * sizeof(struct star_cell *));  
star2 = new star_cell [length+1];
     for (long int i=1; i<=length; i++)
     {
         
         star2[i].cell=star[i].cell;    
         star2[i].mass=star[i].mass; 
         star2[i].sn=star[i].sn; 
         star2[i].sterilized=star[i].sterilized; 
         star2[i].subsection_z=star[i].subsection_z; 
         star2[i].x_coord=star[i].x_coord; 
         star2[i].y_coord=star[i].y_coord; 
         star2[i].z_coord=star[i].z_coord;
         
         star2[i].birth_date=star[i].birth_date;
         star2[i].death_date=star[i].death_date;
         star2[i].ms_lifetime=star[i].ms_lifetime;
         star2[i].sterilization_date=star[i].sterilization_date;
         star2[i].sterilized_count=star[i].sterilized_count;
         star2[i].sterilization_distance=star[i].sterilization_distance;
         
        
      }

  

}


void create_columns()
{

ofstream column_test ("column_test.txt");


column_array = (struct column_properties **)realloc(column_array, (count_xy_glbl+1) * sizeof(struct column_properties *));            
     
   
     for (int l=1; l<=count_xy_glbl; l++)
    {
          column_array[l] = (struct column_properties *)malloc(sizeof(struct column_properties));
          
          column_array[l]->cell_min=1+(l*count_z_glbl)-count_z_glbl;
          column_array[l]->cell_max=l*count_z_glbl;
          
          
          if (l==1)
          {
          column_array[l]->num_stars=
          subsect_param_array[column_array[l]->cell_max]->star_index_max;
          }
          else
          {
         // previous_max=column_array[l-1]->num_stars;   
          column_array[l]->num_stars=
          subsect_param_array[column_array[l]->cell_max]->star_index_max-
          subsect_param_array[column_array[l-1]->cell_max]->star_index_max;
          }
          
          
          
          column_test <<column_array[l]->cell_min<<"," <<column_array[l]->cell_max<<","<<column_array[l]->num_stars<<endl;
    
        
    }    
     
     column_test.close();
    
}

//selects between Salpeter and Kroupa IMFs
float stellar_mass()
{

    return (Salpeter_IMF());  
                
}
      

//produces a stellar mass using the Salpeter IMF
float Salpeter_IMF()
{
      float m_max_fn=m_max;     
      float m_min_fn=m_min;
      float alpha_fn=alpha;
      
      float stellar_mass;
      
      float x = rg.Random();
      
      stellar_mass=
      pow
      (
              (
                  x*
                  (
                      
                      (pow(m_max_fn,(1-alpha_fn)))-
                      (pow(m_min_fn,(1-alpha_fn)))
                      
                  )+
                  pow(m_min_fn,(1-alpha_fn))
              ),
              (1/(1-alpha_fn))
      );
      
      
      
      
      return stellar_mass;

}


long int whole_num_star(float number)
{
    float r_float;
    long int whole_num=(int)number;
    
    float remainder=number-whole_num;
    r_float = rg.Random();
    
    if(r_float<=remainder)
    {
    return (whole_num+1);
    }
    else
    return(whole_num);

}
   


//when you know the CELL index-to the left
//in star1 from the perspective of star2
long int get_star_range_lower_cell_left(int cell)
{
    long int stellar_index;

    //if the cell is the first one in the column, it will divide with a remainder of 1
    if(((cell-count_z_glbl)%count_z_glbl)==1)
    {
    stellar_index=subsect_ptr[cell-count_z_glbl]->star_index_min-(running_total((subsect_ptr[cell]->subsection_xy)-1)); 
    }
    else
    {
    stellar_index=subsect_ptr[cell-count_z_glbl-1]->star_index_min-(running_total((subsect_ptr[cell]->subsection_xy)-1));    
    }
    
    return(stellar_index);
    
}

//when you know the CELL index-to the left
//in star1 from the perspective of star2
long int get_star_range_upper_cell_left(int cell)
{
    long int stellar_index;
    
    //if the cell is the last one in the column, it will divide evenly into the nunber of cells in the column
    if(((cell-count_z_glbl)%count_z_glbl)==0)
    {
    stellar_index=subsect_ptr[cell-count_z_glbl]->star_index_max-(running_total((subsect_ptr[cell]->subsection_xy)-1));
    }
    else
    {
    stellar_index=subsect_ptr[cell-count_z_glbl+1]->star_index_max-(running_total((subsect_ptr[cell]->subsection_xy)-1));    
    }
    
    return(stellar_index);
    
}


//when you know the CELL index-to the right
//which is star2
long int get_star_range_lower_cell_right(int cell)
{
    long int stellar_index;
    int runningtotal;
    
  
    //if the cell is the first one in the column, it will divide evenly
    if(((cell+count_z_glbl)%count_z_glbl)==1)
    {
    stellar_index=subsect_ptr[cell+count_z_glbl]->star_index_min-(running_total(subsect_ptr[cell]->subsection_xy+1));
    }
    else
    {
    stellar_index=subsect_ptr[cell+count_z_glbl-1]->star_index_min-(running_total(subsect_ptr[cell]->subsection_xy+1));    
    }
    
    return(stellar_index);
    
}

//when you know the CELL index-to the right
long int get_star_range_upper_cell_right(int cell)
{
    long int stellar_index;
    
    //if the cell is the last one in the column, it will divide evenly into the nunber of cells in the column
    if(((cell+count_z_glbl)%count_z_glbl)==0)
    {
    stellar_index=subsect_ptr[cell+count_z_glbl]->star_index_max-(running_total(subsect_ptr[cell]->subsection_xy+1));
    }
    else
    {
    stellar_index=subsect_ptr[cell+count_z_glbl+1]->star_index_max-(running_total(subsect_ptr[cell]->subsection_xy+1));    
    }
    
    return(stellar_index);
    
}
//when you know the star index- for the self column
long int get_star_range_lower(int sn_star)
{
      
    long int cell=star1[sn_star].cell;

    
    int lower=0;
    int subsection=star1[sn_star].subsection_z;
    int prev_col_max;

    //if it's the very first star
    if(sn_star==1)
    {
         lower=1;         
     
    }
    //if its the very first column
    else if (subsect_ptr[cell]->subsection_xy==1)
    {
        // the very first cell in the first column
        if (subsection==1)
            {
                 
                 lower=1;              
            }
            //if its not the first cell in the column
            else
            {
                lower=subsect_ptr[cell-1]->star_index_min;   
            }
    }
    //if its not the first column
    //first cell
    else if (subsection==1 && (subsect_ptr[cell]->subsection_xy>1))
    {
         prev_col_max=running_total(subsect_ptr[cell]->subsection_xy);
         lower=subsect_ptr[cell]->star_index_min-prev_col_max;               
    }
    //not the first cell
    else
    {
        prev_col_max=running_total(subsect_ptr[cell]->subsection_xy);
        lower=subsect_ptr[cell-1]->star_index_min-prev_col_max;   
    }
   
   //if the lower value is 0 because the star index has not changed from 0
   //change it to 1
   if (lower==0)
   {
                lower=1;
   }

    
    return(lower);
    
       
}

//when you know the star index- for the self column
long int get_star_range_upper(int sn_star)
{
     
    
    int cell=star1[sn_star].cell;
    int subsection=star1[sn_star].subsection_z;
    long int upper;
    int column=subsect_ptr[cell]->subsection_xy;
    
    
    //if its the first column
    if (subsect_ptr[cell]->subsection_xy==1)
    {

        if (subsection==count_z_glbl)
            {
                 
                 upper=subsect_ptr[cell]->star_index_max;            
            }
            else
            {
                upper=subsect_ptr[cell+1]->star_index_max;   
            }                                
    
    }
    
    //if its the last cell- not in the first column
    else if (subsection==count_z_glbl)
    {
         upper=subsect_ptr[cell]->star_index_max-(running_total(column));               
    }
    else
    {
         upper=subsect_ptr[cell+1]->star_index_max-(running_total(column)); 
    }
   
    
    return(upper);
    
       
}

int SN_star(float mass)
{
int sn=0;
float typeII_mass_fn=typeII_mass;

     if (mass>=typeII_mass_fn)
     {
         sn=1;
     }        

     return (sn);
      
}

float degree_generator(float deg_min, float deg_max)
{
      //int r_int=0;
      int32 r_int=0;
      float r_float=0;

    
     //produces a float between 0 and 1
     r_float = rg.Random();
     r_int = rg.IRandomX(deg_min,deg_max-1);
     
     
     return (r_float+r_int);
      
}


float random_star_generator_Z(float z_min, float z_max)
{
      
      int r_int;
      float r_float;
    
     //produces a float between 0 and 1
     r_float = rg.Random();
     r_int = rg.IRandom(z_min,z_max-1);
     
     return (r_float+r_int);
      
}


float random_star_generator_R (float deg_min, float deg_max, float r_min, float r_max)
{
      float R=0;
      float w;
      //int r_int;
      float r_float;
      
     //produces a float between 0 and 1
     r_float = rg.Random();
     w=r_float;
     
     
     R=sqrt(
     w*((r_max*r_max)-(r_min*r_min))
     +(r_min*r_min)
     );
     
     //cout <<"\n r_float: "<<r_float; 
     //cout << "\n R: "<<R;
     
      return R;     
}

void create_subsections()
{
     cout <<"\nCreating Model Structure (Cells, Columns etc.)";
    int count_z=1;
    int count_xy=1;
    int number_zs_iterated=0;
    int number_xys_interated=0;
    
    
    //struct subsection_parameters **subsect_param_array = NULL;
    
    
    float mw_radius_fn=mw_radius;
    float mw_height_fn=mw_height;
    float mw_angle_fn=mw_angle;
    float subsection_slice_radial_distance_fn=subsection_slice_radial_distance; 
    float subsection_slice_h_z_distance_fn=subsection_slice_h_z_distance;
    float pi_fn=PI;
    float area;
   
 

    int xy_subsections= ceil(mw_radius_fn/subsection_slice_radial_distance_fn);
    int z_subsections= ceil(mw_height_fn/subsection_slice_h_z_distance_fn);
    
    //height of mw in z-direction after ceiling function performed on the number of slices
    float mw_height_real=z_subsections*subsection_slice_h_z_distance_fn;
    
    float p_min=0;
    float p_max=subsection_slice_radial_distance_fn;
    
    //height above the mid-plane
    float mw_height_above_mp=mw_height_real/2;

   int number_density_type_fn=number_density_type;
   
   float h_R_fn= h_R; 
   float Z_thin_fn=Z_thin; 
   float Z_thick_fn= Z_thick; 
   float n0_fn=n0; 
   


   float total_volume=0;
  
  //allocate memory for each subsection
  //assign each subsection a z-identifier and xy-plane identifier
   while(count_xy!=xy_subsections+1)
   {
           
           p_min=number_xys_interated*subsection_slice_radial_distance_fn;
           p_max=subsection_slice_radial_distance_fn+(number_xys_interated*subsection_slice_radial_distance_fn);
           
          
           
           mw_height_above_mp=(z_subsections*subsection_slice_h_z_distance_fn)/2;                         
    
           
           //sets subsection parameters for the z-subsection
           while(count_z!=z_subsections+1)
           {
          
            subsect_param_array = (struct subsection_parameters **)realloc(subsect_param_array, (count_z+(number_zs_iterated*z_subsections)+ 1) * sizeof(struct subsection_parameters *));
            //allocate memory for one subsection_parameter structure
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]= (struct subsection_parameters*)malloc(sizeof(struct subsection_parameters));
           
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->subsection_xy=count_xy;
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->subsection_z=count_z;
            
           
            //set the minimum and maximum angle
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->angle_min=0;
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->angle_max=mw_angle;
           
           //set the max and min heights above/below the midplane
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_max=mw_height_above_mp;
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_min=mw_height_above_mp-subsection_slice_h_z_distance_fn;
            
            mw_height_above_mp=mw_height_above_mp-subsection_slice_h_z_distance_fn;
            
            
            //set the minimum and maximum points for the cylindrical coordinates
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_min=p_min;
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_max=p_max;
            
            
            
            //set the midpoint to this subsection in xy and z in cartesian coordinates
            
            
            //xy-radial distance
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->midpoint_xy_plane_dist=
            (
            (
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_max
            -subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_min
            )
            /2)
            +
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_min
            ;
            
            //z- make the value positive (height above midplane) (if it's a negative) for the stellar density formula
            //one exception is when one value is +ve and other -ve, (when the z-midpoint is 0-at the midplane)
            if ((subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_max+subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_min)==0)
            {
               subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->mp_z_coord_height_above_midplane=0;                                                                                                                                                      
            }
            else
            {
                subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->mp_z_coord_height_above_midplane=
                (
                (
                sqrt((subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_max)*(subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_max))+
                sqrt((subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_min)*(subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->z_min))
                )
                /
                (2.0)
                );
            }
            //determine the volume of the subsection
            //its the volume of the p_max-p_min (overlapping circles)
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->volume=
            (
            (mw_angle_fn/360)*(pi_fn)*
            (subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_max)*
            (subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_max)*
            (subsection_slice_h_z_distance_fn)
            )
           
            -
            (
            (mw_angle_fn/360)*(pi_fn)*
            (subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_min)*
            (subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->radial_distance_min)*
            (subsection_slice_h_z_distance_fn)
            )
            ;          
            
            total_volume=total_volume+subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->volume;
            //determine the number of stars in the subsection now that the volume is known
            //use stellar density in carroll and ostlie or Dehnon & Binney
               
            //number density carroll and ostlie    
            if (number_density_type_fn==1)
            {   
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->num_stars=whole_num_star(
            (subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->volume)*
            (
            (n0_fn)*(
            exp(-subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->mp_z_coord_height_above_midplane/Z_thin_fn)+ 
            0.085*exp(-subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->mp_z_coord_height_above_midplane/Z_thick_fn) 
            
            )*
            exp(-subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->midpoint_xy_plane_dist/h_R_fn)
            ))
            ;
            }
            
           
                    
                    
  
            
            if((subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->num_stars+num_stars_glbl)>=(num_stars_glbl+1))
            {
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->star_index_min=num_stars_glbl+1; 
            }
            else
            {
                subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->star_index_min=num_stars_glbl; 
            }
            
            num_stars_glbl=num_stars_glbl+subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->num_stars; 
            
            subsect_param_array[count_z+(number_zs_iterated*z_subsections)]->star_index_max=num_stars_glbl; 
            
            count_z++;
            
          }    
          
          count_z=1;
          count_xy++;
          number_zs_iterated++;
          number_xys_interated++;       
   }
   
   cout <<"\n num stars: "<<setprecision(15)<<num_stars_glbl;
  
  //global pointer to structure of subsection parameters
  subsect_ptr=subsect_param_array;
  
  //setting the number of z and xy subsections
  count_z_glbl=z_subsections;
  count_xy_glbl=xy_subsections;
  
} //end of function




//gathers stats after each column is processed- radial
void print_stats()
{
    
     int max_len_ster_hist_fn=max_len_ster_hist;
     long int total_sterilizations=0;
     long long typeII_sterilizations=0;
     long long typeIa_sterilizations=0;
  
     float radial_dist_cell_size=subsection_slice_radial_distance;
     stats_counter++;
     
     int planet_habitable_for_4Gyr=0;

                       
     for (int i=1; i<=column_array[stats_counter]->num_stars; i++)
     {
     
         //make sure that a planet existing around a star is 0 to start-resolve later
       planet_habitable_for_4Gyr=0;
        
         if ((star1[i].sterilized)==1 || (star1[i].sterilized)==2)
         {
              total_sterilizations++;    //star was at least sterilized once                                          
         }
         
 
     } //end of for loop that iterates through every column

     stats[stats_counter] = (struct statistics *)malloc(sizeof(struct statistics));
     stats[stats_counter]->total_sterilizations=total_sterilizations; 
     stats[stats_counter]->distance=subsect_ptr[(column_array[stats_counter]->cell_min)]->midpoint_xy_plane_dist;
     stats[stats_counter]->total_stars=column_array[stats_counter]->num_stars;
     
} //end of function

//outputs stats to file
void output_stats()
{
     ofstream GHZ_output("GHZ_output.txt", ios::out);
     GHZ_output << "column_number, radial distance, total stars, total sterilizations- at least sterilized once\n";
     
  
     for (int i=1; i<=count_xy_glbl; i++)
     {
          GHZ_output << i <<","
          <<stats[i]->distance<<","
          <<stats[i]->total_stars*360<<","
          <<stats[i]->total_sterilizations*360
          <<endl;
                          
     }

    GHZ_output << "The following results are multiplied by 360 such that they apply to the entire Galaxy\n"<< "Number of TypeII: "<< num_sn_glbl*360<<" Number of TypeIa: "<<num_sn_Ia_glbl*360<<" Total Sterilizations SNII: "<<total_sterilizations_II*360<<" Total Sterilizations Ia: "<<total_sterilizations_Ia*360<<" Total stellar mass: "<<total_stellar_mass_glbl*360;
    GHZ_output <<"Number of SNII: "<< num_SNII_specific*360;
    GHZ_output <<"\nNumber of SNIa: "<< num_SNIa_specific*360;
    
float slice_radial_distance=subsection_slice_radial_distance;
float slice_h_z_distance=subsection_slice_h_z_distance; 
int z_cell_processed=0;
     
     
      GHZ_output.close();
}
