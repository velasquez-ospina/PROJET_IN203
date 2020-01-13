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
    const int nb_ants = 2*dims.first*dims.second; // Nombre de fourmis
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
    // On va créer toutes les fourmis dans le nid :
    std::vector<ant> ants;
    ants.reserve(nb_ants);
    size_t food_quantity = 0;
    start2 = std::chrono::system_clock::now();
    for ( size_t i = 0; i < nb_ants; ++i )
        ants.emplace_back(pos_nest, life);
    // On crée toutes les fourmis dans la fourmilière.
    pheronome phen(laby.dimensions(), pos_food, pos_nest, alpha, beta);
    const unsigned int buffer_size = dims.first*dims.second+1+2*nb_ants;

	MPI_Init( &nargs, &argv );                 
	MPI_Comm_size(MPI_COMM_WORLD, &nbp);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    if (rank == 0){

        double buffer[buffer_size];
        gui::context graphic_context(nargs, argv);
        gui::window& win =  graphic_context.new_window(h_scal*laby.dimensions().second,h_scal*laby.dimensions().first+266);
        display_t displayer( laby, phen, pos_nest, pos_food, ants, win );

        gui::event_manager manager;
        manager.on_key_event(int('q'), [] (int code) { exit(0); });
        manager.on_key_event(int('t'), [] (int code) { std::cout << "Temps pour le loop advance_time: " << time1.count()  << std::endl;});
        manager.on_display([&] { displayer.display(food_quantity); win.blit(); });
        manager.on_idle([&] () { 
            displayer.display(food_quantity);
            win.blit(); 
            if (food_quantity > 100 && Switch_end){

                end2 = std::chrono::system_clock::now();
                std::chrono::duration<double> time2 = end2-start2;
                std::cout << "temps pour 100 de nourriture : " << time2.count()  << std::endl;
                Switch_end = false;
            }
            MPI_Recv(buffer,buffer_size,MPI_DOUBLE,MPI_ANY_TAG,101,MPI_COMM_WORLD, &status);

            food_quantity = buffer[0];
            //for (unsigned int i= 0; i < dims.first;i++)
              //  for (unsigned int j = 0; j < dims.second ; j++)

            for (int k = 0,l = 0; k < nb_ants; k++,l+=2)
            {
                position_t ant_pos;
                ant_pos.first =buffer[dims.first*dims.second+1+l];
                ant_pos.second =buffer[dims.first*dims.second+2+l];
                ants[k].set_position(ant_pos);
            } 
        });
        manager.loop();
    }else{
        while(1){
            advance_time(laby, phen, pos_nest, pos_food, ants, food_quantity);
            double buffer[buffer_size];
            buffer[0]=food_quantity;
            for (unsigned int i= 0; i < dims.first;i++)
                for (unsigned int j = 0; j < dims.second ; j++)
                    buffer[dims.first*i+j+1] = phen(i,j);

            for (int k = 0,l = 0; k < nb_ants; k++,l+=2)
            {
                position_t ant_pos=ants[k].get_position();
                buffer[dims.first*dims.second+1+l]= ant_pos.first;
                buffer[dims.first*dims.second+2+l]= ant_pos.second;
            }

            MPI_Send(buffer,buffer_size,MPI_DOUBLE,0,101,MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    return 0;
}