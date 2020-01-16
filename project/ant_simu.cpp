#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include "labyrinthe.hpp"
#include "ant.hpp"
#include "pheronome.hpp"
# include "gui/context.hpp"
# include "gui/colors.hpp"
# include "gui/point.hpp"
# include "gui/segment.hpp"
# include "gui/triangle.hpp"
# include "gui/quad.hpp"
# include "gui/event_manager.hpp"
# include "display.hpp"
# include <mpi.h>

/*
std::chrono::time_point<std::chrono::system_clock> start1, end;
  start1 = std::chrono::system_clock::now();

  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start1;
  std::cout << "Temps assemblage vecteurs : " << elapsed_seconds.count()  << std::endl;


*/

    std::chrono::time_point<std::chrono::system_clock> start1, start2, start3, start4, end1, end2, end3, end4;
    std::chrono::duration<double> time1, time2, time3, time4;
void advance_time( const labyrinthe& land, pheronome& phen, 
                   const position_t& pos_nest, const position_t& pos_food,
                   std::vector<ant>& ants, std::size_t& cpteur)
{
    start1 = std::chrono::system_clock::now();
    for ( size_t i = 0; i < ants.size(); ++i ){
        ants[i].advance(phen, land, pos_food, pos_nest, cpteur);}
    end1 = std::chrono::system_clock::now();
    time1 = end1-start1;
    phen.do_evaporation();
    phen.update();
}


int main(int nargs, char* argv[])
{
    int nbp,rank;
    bool Switch_end = true;
    const dimension_t dims{32,64};// Dimension du labyrinthe
    const std::size_t life = int(dims.first*dims.second);
    const double eps = 0.75;  // Coefficient d'exploration
    const double alpha=0.97; // Coefficient de chaos
    //const double beta=0.9999; // Coefficient d'évaporation
    const double beta=0.999; // Coefficient d'évaporation
    labyrinthe laby(dims);

    // Location du nid
    position_t pos_nest{dims.first/2,dims.second/2};
    // Location de la nourriture
    position_t pos_food{dims.first-1,dims.second-1};
    // Définition du coefficient d'exploration de toutes les fourmis.
    ant::set_exploration_coef(eps);
    MPI_Init( &nargs, &argv );                 
	MPI_Comm_size(MPI_COMM_WORLD, &nbp);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    // On va créer toutes les fourmis dans le nid :
    std::vector<ant> ants;
    std::vector<ant> local_ants;
    size_t food_quantity = 0;
    start2 = std::chrono::system_clock::now();

    const unsigned int nb_ants = 2*dims.first*dims.second; // Nombre de fourmis
    const unsigned int local_nb_ants = nb_ants/(nbp-1);
    const int buffer_size = dims.first*dims.second+1+2*nb_ants;
    const int buffer_ants_size = 2*local_nb_ants;
    const int buffer_phen_size = dims.first*dims.second;
    double buffer[buffer_size]; //to communicate between process 0 and 1
    
    ants.reserve(nb_ants);
    local_ants.reserve(local_nb_ants);
    pheronome phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);

    MPI_Group world_group, group_0_1, group_1_n;
    MPI_Comm comm_0_1, comm_1_n;

    int ranks_0_1[2] = {0, 1};
    const unsigned int process = nbp-1;
    int ranks_1_n[process];
   
    for(int i = 1; i < nbp; i++){
        ranks_1_n[i - 1] = i; 
    }

    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    if(rank == 1 || rank == 0){
        MPI_Group_incl(world_group, 2, ranks_0_1, &group_0_1);
        MPI_Comm_create_group(MPI_COMM_WORLD, group_0_1, 0, &comm_0_1);
    }
    if(rank >= 1){
        MPI_Group_incl(world_group, nbp - 1, ranks_1_n, &group_1_n);
        MPI_Comm_create_group(MPI_COMM_WORLD, group_1_n, 0, &comm_1_n);   
    }
    for ( unsigned int i = 0; i < nb_ants; ++i ){
        ants.emplace_back(pos_nest, life);
    }


    if (rank == 0){

        std::vector<double> phen_buffer;
        gui::context graphic_context(nargs, argv);
        gui::window& win =  graphic_context.new_window(h_scal*laby.dimensions().second,h_scal*laby.dimensions().first+266);
        display_t displayer( laby, phen, pos_nest, pos_food, ants, win );

        gui::event_manager manager;
        manager.on_key_event(int('q'), [] (int code) { exit(0); });
        manager.on_key_event(int('f'), [&food_quantity] (int code) { std::cout << "Food : " << food_quantity << std::endl;});
        manager.on_key_event(int('t'), [] (int code) { std::cout << "Temps pour le loop advance_time: " << time1.count()  << std::endl;});
        manager.on_display([&] { displayer.display(food_quantity); win.blit(); });
        manager.on_idle([&] () { 
            displayer.display(food_quantity);
            win.blit(); 
            if (food_quantity > 10000 && Switch_end){

                end2 = std::chrono::system_clock::now();
                std::chrono::duration<double> time2 = end2-start2;
                std::cout << "temps pour 10000 de nourriture: " << time2.count()  << std::endl;
                Switch_end = false;
            }
                MPI_Recv(buffer,buffer_size,MPI_DOUBLE,1,101,comm_0_1, &status);

                food_quantity = buffer[0];
                for (unsigned int i= 0; i < dims.first;i++){
                    for (unsigned int j = 0; j < dims.second ; j++){
                        phen_buffer.emplace_back(buffer[i*dims.second+j+1]);
                    }
                }
                phen.Swap_copy(phen_buffer);
                phen_buffer.resize(0);

                for (unsigned int k = 0,l = 0; k < nb_ants; k++,l+=2) {
                    position_t ant_pos;
                    ant_pos.first =buffer[buffer_phen_size+1+l];
                    ant_pos.second =buffer[buffer_phen_size+2+l];
                    ants[k].set_position(ant_pos);
                } 
            
        });
        manager.loop();
    }else if (rank >= 1){
        double buffer_phen[buffer_phen_size];
        double buffer_ants[buffer_ants_size];
        int buffer_food;

        //pour recevoir les fourmis
        double recv_ants[buffer_ants_size*(nbp-1)];

        //pour recevoir les phens
        double recv_phens[buffer_phen_size];
        std::vector<double> phen_buffer;
        for ( unsigned int i = 0; i < local_nb_ants; ++i )
            local_ants.emplace_back(pos_nest, life);
        // On crée toutes les fourmis dans la fourmilière.
        while(1){
            advance_time(laby, phen, pos_nest, pos_food, local_ants, food_quantity);


            MPI_Reduce (&food_quantity, &buffer_food, 1, MPI_INT, MPI_SUM,0, comm_1_n );
            //packing phens
            for (unsigned int i= 0; i < dims.first;i++){
                for (unsigned int j = 0; j < dims.second ; j++){
                    buffer_phen[dims.second*i+j] = phen(i,j);
                }
            }
            MPI_Allreduce (buffer_phen, recv_phens, buffer_phen_size, MPI_DOUBLE,MPI_MAX, comm_1_n );
            //Unpacking phens
            for (unsigned int i= 0; i < dims.first;i++){
                for (unsigned int j = 0; j < dims.second ; j++){
                    phen_buffer.emplace_back(recv_phens[i*dims.second+j]);
                }
            }        
            phen.Swap_copy(phen_buffer);
            phen_buffer.resize(0);

            //packing ants
            for (unsigned int k = 0,l = 0; k < local_nb_ants; k++,l+=2)
            {
                position_t ant_pos=local_ants[k].get_position();
                buffer_ants[l]= ant_pos.first;
                buffer_ants[l+1]= ant_pos.second;
            }
            MPI_Gather(buffer_ants, buffer_ants_size, MPI_DOUBLE, recv_ants, buffer_ants_size, MPI_DOUBLE, 0, comm_1_n);

            if(rank == 1){
                buffer[0]=(double)buffer_food;    //first position of the buffer is the food quantity
                //packing the new phens in the buffer that goes to the display
                for (unsigned int i= 0; i < dims.first;i++){
                    for (unsigned int j = 0; j < dims.second ; j++){
                        buffer[dims.second*i+j+1] = phen(i,j);
                    }
                }
                 //packing the ants in the buffer
                for (unsigned int k = 0,l = 0; k < nb_ants; k++,l+=2)
                {
                    buffer[buffer_phen_size+1+l]=recv_ants[l];
                    buffer[buffer_phen_size+2+l] =recv_ants[l+1];
                } 
                MPI_Send(buffer,buffer_size,MPI_DOUBLE,0,101,comm_0_1);
            }

        }

    }
    MPI_Finalize();
    return 0;
}