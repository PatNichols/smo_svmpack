#include "putils_fs.h"
#include "putils_program_options.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>

void create_command_line(int& argc,char ***argv)
{
    std::string mess[7] =
    {
        "-option1=value1",
        "-option3",
        " value 3",
        "--option11=value11",
        "--option5",
        "value5",
        "option12"
    };
    argc = 7; 
    *argv = new char*[argc];
    for (auto i=0;i<(argc);++i) {
        (*argv)[i] = Strdup(mess[i].c_str());
    }
}

void set_config()
{
    std::ofstream is("test.config");
    is << "option8 value8\n";
    is << "option9\n";
    is << "option6 value6\n";
    is.close();
}

void set_env()
{
    setenv("TEST_VALUE7","value7",1);
}


void set_options()
{
    const char * name = "test.options";
    std::ofstream out(name);
    out << "options1 \"options 1\" required" << "\n";
    out << "options2 \"options 1\" value2 required" << "\n";
    out << "options3 \"options 1\" value3" << "\n";
    out << "options4 required " << "\n";
    out << "options5 value5 required" << "\n";
    out << "options6 value6"<< "\n";
    out << "options7 value7" << "\n";
    out << "options8" << "\n";
    out << "options7" << "\n";
    out << "options8" << "\n";
    out << "options8" << "\n";
    out << "options10" << "\n";
    out << "options11" << "\n";
    out << "options12" << "\n";
    out.close();
}


int main()
{
  try {
    char **argv;
    int argc;
    set_env();
    set_options();
    set_config();
    create_command_line(argc,&argv);
    putils::ProgramOptions popts;
    popts.addOptionsFromFile("test.options");
    std::cout << popts << "\n";
    std::cout <<" ------------------------------------\n";
    popts.parseCommandLine(argc,argv);
    std::cout << popts << "\n";
    std::cout <<" ------------------------------------\n";
    popts.parseFile("test.config");
    std::cout << popts << "\n";
    std::cout <<" ------------------------------------\n";
    popts.parseEnv("TEST");
    std::cout << popts << "\n";
    std::cout <<" ------------------------------------\n";    
  } catch (std::exception& e) {
     e.what();
     exit(-1);  
  }
}