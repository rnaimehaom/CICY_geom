//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Main function (program start)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include <string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "App_Network_3D.h"
#include "App_Content_Dist.h"
#include "App_Network_From_Abaqus.h"

int main(int argc, char** argv)
{
    //Read input file name into in_file
    string in_file;
    if(argc > 1)    in_file = argv[1];
    else
    {
        cout << "The input file name is:  ";
        in_file = "input.dat";
        cout << in_file << endl;
    }
    //Open the input file
    ifstream infile;
    infile.open(in_file.c_str());
    if(!infile) { hout << "Failed to open input file: "  << in_file << endl;  return 0; }

    //Read output file name into out_file
    string out_file;
    if(argc > 2)    out_file = argv[2];
    else
    {
        cout << "The output file name is:  ";
        out_file = "output.dat";
        cout << out_file << endl;
    };
    //Open the output stream
    if(out_file.size()>0) open_deffo_stream( (char*)out_file.c_str() );

    //----------------------------------------------------------------------
    //Identification Tag

    hout<<endl;
    hout<<"******************************************"<<endl;
    hout<<"            3D Geometric model      "<<endl;
    hout<<"                                          "<<endl;
    hout<<"  Propriety of CICY/Unidad de Materiales    "<<endl;
    hout<<"              faviles@cicy.mx       "<<endl;
    hout<<"******************************************"<<endl;
    hout<<endl;
    hout<<endl;

    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //Call for application cases

    //Time markers for total simulation
    time_t it_begin, it_end;
    it_begin = time(NULL);

    //----------------------------------------------------------------------
    //Input file reader
    hout<<"======================================================" << endl;
    hout<<"Reading input file......"<<endl;
    Input *Init = new Input;
    //Initialize input data with default values
    if(Init->Data_Initialization())
    {
        //Read the input file
        if (!Init->Read_input_file(infile)) {
            hout << "Error when reading input file." << endl;
            return 0;
        }
    }
    else {
        hout << "Error when initializing input data." << endl;
        return 0;
    }
    it_end= time(NULL);
    hout<<"Input file read in "<<(int)(it_end-it_begin)<<" secs."<<endl;

    //Close input file
    infile.close();

    //----------------------------------------------------------------------
    //Check which application is being called
    if(Init->app_name.str=="3D_Electrical_Network")
    {
        //----------------------------------------------------------------------
        //Application to create a 3D network of nanotubes
        App_Network_3D *Network3D = new App_Network_3D;
        int count = Init->simu_para.sample_num;
        //Implement all samples
        for(int i=1; i<=count; i++) {
            if(!Network3D->Generate_nanoparticle_resistor_network(Init)) {
                hout<<"Error in Generate_nanoparticle_resistor_network in sample "<<i<<endl;
                return 0;
            }
        }
            
        delete Network3D;
    }
    else if(Init->app_name.str=="Content_Distribution") {
        
        //Create an object to call the application that calculates the content
        //at each observation window
        App_Content_Dist *ContentDist = new App_Content_Dist;
        if (!ContentDist->Calculate_content_on_each_window(Init)) {
            hout<<"Error in Calculate_content_on_each_window"<<endl;
            return 0;
        }
        delete ContentDist;
    }
    else if (Init->app_name.str == "Network_From_Abaqus")
    {
        //Create an object to call the application that reads a network from an Abaqus database
        App_Network_From_Abaqus* NetworkAbqs = new App_Network_From_Abaqus;
        if (!NetworkAbqs->Nanoparticle_resistor_network_from_odb(Init)) {
            hout << "Error in Nanoparticle_resistor_network_from_odb" << endl;
            return 0;
        }
        delete NetworkAbqs;
    }
    else {
        //Delete Init object and terminate with error
        hout<<"Error: Invalid application name. Input was:"<<Init->app_name.str<<endl;
        delete Init;
        return 0;
    }
    
    //----------------------------------------------------------------------
    //Delete the object with the input data
    delete Init;

    //----------------------------------------------------------------------
    //Time markers for total simulation
    it_end = time(NULL);
    hout<<endl;
    hout<<"*******************************************************************"<<endl;
    hout<<"    End of simulation "<<endl;
    hout<<"    Simulation took "<<(int)(it_end-it_begin)<<" secs."<<endl;
    hout<<"*******************************************************************"<<endl;
    

    //----------------------------------------------------------------------
    //Close the output stream
    close_deffo_stream();
    return 1;

}
//===========================================================================
