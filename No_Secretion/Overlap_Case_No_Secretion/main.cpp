/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	// load and parse settings file(s)
	
	bool XML_status = false; 
	char copy_command [1024]; 
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }
	
	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	
	setup_microenvironment(); // modify this in the custom code 
	
    Microenvironment coarse_well;
    coarse_well.name = "coarse_well";
    coarse_well.spatial_units = "micron";
    coarse_well.mesh.units = "micron";
    coarse_well.time_units = "min";
    
    coarse_well.set_density( 0 , "oxygen", "mM", 108000.0 , 0.00 ); //108000
    coarse_well.add_density( "glucose", "mM", 30000.0 , 0.0 ); //30000
    coarse_well.add_density( "chemokine", "mM", 100000 , 0.0);
    coarse_well.resize_space( 100, 1 , 1 );
    
    double dx = 32;
    double dy = 2880;
    double dz = 2880;
    
    coarse_well.resize_space( 224.0, 5120.0 , -dy/2.0+16 , dy/2.0+16 , -dz/2.0+16 , dz/2.0+16 , dx, dy, dz );
    std::vector<double> dirichlet_condition = { 0 , 0, 0, 0 };

    coarse_well.set_substrate_dirichlet_activation(0,false);
    coarse_well.set_substrate_dirichlet_activation(1,false);
    coarse_well.set_substrate_dirichlet_activation(2,false);
    
    coarse_well.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;   
    
    for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
    {
        coarse_well(m)[0]= 0.285; // oxygen  //0.285
        coarse_well(m)[1]= 16.897255; // glucose //16.897255
        coarse_well(m)[2]= 0; //chemokine
    }
	
	
	for ( int m = 0; m < coarse_well.mesh.voxels.size() ; m++)
    {
		double coarMic_y = coarse_well.mesh.voxels[m].center[0];
		//std::cout << coarMic_y << std::endl;
		if (coarMic_y == 240)
		{
			std::cout << "changing voxels with y = 240" << std::endl;
			coarse_well(m)[0]= 0 ; // oxygen
			coarse_well(m)[1]= 0; // glucose
			coarse_well(m)[2]= 0; //chemokine
		}
		
	}
    
    coarse_well.display_information( std::cout );
    coarse_well.write_to_matlab("output/output00000000_microenvironment1.mat");
    
    /*
    Microenvironment transfer_region;
    transfer_region.name = "transfer_region";
    transfer_region.spatial_units = "micron";
    transfer_region.mesh.units = "micron";
    transfer_region.time_units = "min";
    
    transfer_region.set_density( 0 , "oxygen", "mmHg", 108000 , 0.00 );
    transfer_region.add_density( "glucose", "mM", 30000 , 0.0 );
    transfer_region.add_density( "chemokine", "mM", 100000 , 0.0);
    transfer_region.resize_space( 100, 1 , 1 );
    transfer_region.resize_space( 224.0, 288.0 , -dy/2.0+16 , dy/2.0+16 , -dz/2.0+16 , dz/2.0+16 , dx, dy, dz );
    
    transfer_region.diffusion_decay_solver = diffusion_decay_solver__constant_coefficients_LOD_1D;   
    
    for ( int m = 0; m < transfer_region.mesh.voxels.size() ; m++)
    {
        transfer_region(m)[0]=38; // oxygen
        transfer_region(m)[1]=16.897255; // glucose
        transfer_region(m)[2]=0; //chemokine
    }
    
    transfer_region.display_information( std::cout );
    transfer_region.write_to_matlab("output/output00000000_microenvironment2.mat");    
    */
    
	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 32; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/* Users typically start modifying here. START USERMODS */ 
	
	create_cell_types();
	
	setup_tissue();

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function; 
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 
	
	display_citations(); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	// main loop 
    
    std::vector<double> v = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 , 20.0 };
	int i = 0;
    int j = 0;
	try 
	{		
        
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1 * diffusion_dt )
		{
            //std::cout << "NEXT TIME WILL BE SAVED AT : " << v[i] << std::endl;
			// save data if it's time. 
			if( fabs( PhysiCell_globals.current_time - v[i] ) < 0.01 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
                    
                    sprintf( filename , "%s/output%08u_microenvironment1.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
                    coarse_well.write_to_matlab(filename);
                    
                    //sprintf( filename , "%s/output%08u_microenvironment2.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
                    //transfer_region.write_to_matlab(filename);
				}
				
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
                i += 1; 
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - v[j]  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
                    j += 1;
				}
			}

			sprintf( filename , "%s/output%08u_microenvironment0_before_diffusion.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			microenvironment.write_to_matlab(filename);
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
            
			sprintf( filename , "%s/output%08u_microenvironment0_after_diffusion.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			microenvironment.write_to_matlab(filename);
            
            // coarsening (transfer_region right hand_side)
            std::vector<double> v = {0, 0, 0};
            for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
            {
                double mic_cen_y = microenvironment.mesh.voxels[m].center[1];
                if (mic_cen_y == 240)
                { 
                    v[0]+=microenvironment(m)[0]*microenvironment.mesh.voxels[m].volume; //oxygen
                    v[1]+=microenvironment(m)[1]*microenvironment.mesh.voxels[m].volume; //glucose
                    v[2]+=microenvironment(m)[2]*microenvironment.mesh.voxels[m].volume; //chemokine
                }
            }
            v[0]=v[0]/coarse_well.mesh.voxels[0].volume;
            v[1]=v[1]/coarse_well.mesh.voxels[0].volume;
            v[2]=v[2]/coarse_well.mesh.voxels[0].volume;

			std::vector<double> fine_overlap = {v[0], v[1], v[2]};

			sprintf( filename , "%s/output%08u_microenvironment1_before_coarsening.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			coarse_well.write_to_matlab(filename);
			
			
			coarse_well(0)[0] = v[0];
            coarse_well(0)[1] = v[1];
            coarse_well(0)[2] = v[2];
			
			sprintf( filename , "%s/output%08u_microenvironment1_after_coarsening.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			coarse_well.write_to_matlab(filename);
            
			std::vector<double> coarse_overlap_before_diffusion = {coarse_well(0)[0], coarse_well(0)[1], coarse_well(0)[2]};

            coarse_well.simulate_diffusion_decay(diffusion_dt);
			
			sprintf( filename , "%s/output%08u_microenvironment1_after_diffusion.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			coarse_well.write_to_matlab(filename);
            
			int coarse_well_voxel_number = coarse_well.mesh.voxels.size();
            // Dirichlet Boundary Condition
            coarse_well(coarse_well_voxel_number-1)[0] = 0.285;
            int y_240 = 0;
			std::vector<double> coarse_overlap_after_diffusion = {coarse_well(0)[0], coarse_well(0)[1], coarse_well(0)[2]};

			sprintf( filename , "%s/output%08u_microenvironment1_after_diffusion_with_DC.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			coarse_well.write_to_matlab(filename);

            // right side overwrite
            //std::cout << y_240 << std::endl;
            double oxy_diff = coarse_overlap_after_diffusion[0];// - coarse_overlap_before_diffusion[0];
            double glu_diff = coarse_overlap_after_diffusion[1];// - coarse_overlap_before_diffusion[1];
            double chem_diff = coarse_overlap_after_diffusion[2];// - coarse_overlap_before_diffusion[2];
			
            sprintf( filename , "%s/output%08u_microenvironment0_before_projection.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			microenvironment.write_to_matlab(filename);
            
            for ( int m = 0; m < microenvironment.mesh.voxels.size() ; m++)
            {
                double mic_cen_y = microenvironment.mesh.voxels[m].center[1];
                if (mic_cen_y == 240)
                { 
                    //std::cout << "Glucose difference per voxel  : "  <<glu_diff << std::endl; 
                    microenvironment(m)[0] = oxy_diff/8100; //oxygen
                    microenvironment(m)[1] = glu_diff/8100; //glucose
                    microenvironment(m)[2] = chem_diff/8100; //chemokine
                }

            }
            
			
			sprintf( filename , "%s/output%08u_microenvironment0_after_projection.mat" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );      
			microenvironment.write_to_matlab(filename);
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			PhysiCell_globals.current_time += diffusion_dt;
            
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 


	return 0; 
}
