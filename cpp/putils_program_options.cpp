#include "putils_program_options.hpp"
namespace putils
{
void ProgramOptions::addOptionsFromFile( const char * name)
{	
  try { 
    std::ifstream in(name);
    if (!in) throw putils::OpenFileError(__FUNCTION__,name,"r");
    std::string sline;
    std::string delims(" \n\0");
    std::vector<std::string> tokens;
    std::string desc;
    while (in) {
        desc.clear();
        if (!std::getline(in,sline)) {
            if ( in.eof() ) break;
            throw putils::ReadError(__FUNCTION__,name);
        }
        if (sline.size() == 0) {
            throw putils::FormatError(__FUNCTION__,name," blank line");
        }
        if (sline[0] == '#') {
            continue;
        }
        size_t p = sline.find_first_of('\"');
        if ( p != std::string::npos) {
            size_t q = sline.find_last_of('\"');
            if ( p == q ) {
                throw putils::FormatError(__FUNCTION__,name,"unterminated quoted string");
            }
            desc = sline.substr(p+1,q-p-1);
            std::string right = sline.substr(0,p);
            std::string left = sline.substr(q+1);
            sline = right;
            sline += left;
        }
        size_t ntokens = putils::tokenize_string(sline,delims,tokens);
        switch (ntokens)
        {
        case 0:
            throw putils::FormatError(__FUNCTION__,name,"not enough tokens on line");
        case 1:
            addOption(tokens[0],desc,std::string(),0);
            break;
        case 2:
            if ( tokens[1].compare("required") == 0) {
                addOption(tokens[0],desc,std::string(),1);
            }else{
                addOption(tokens[0],desc,tokens[1],0);
            }
            break;
        case 3:
            if ( tokens[2].compare("required") == 0) {
                addOption(tokens[0],desc,tokens[1],1);
            }else{
                throw putils::FormatError(__FUNCTION__,name,"too many tokens on line");
            }
            break;
        default:
            throw putils::FormatError(__FUNCTION__,name,"too many tokens on line");
        }
    }
    in.close();
  } catch (std::exception& e) {
     e.what();
     putils_error_function();
  }
}

void ProgramOptions::parseFile( const char * name)
{	
    std::ifstream in(name);
    if (!in) throw putils::OpenFileError(__FUNCTION__,name,"r");
    std::string sline;
    std::string delims(" =\n");
    std::vector<std::string> tokens;
    while (in)
    {
        if (!std::getline(in,sline)) {
            if ( in.eof() ) break;
            throw putils::ReadError(__FUNCTION__,name);
        }
        if (sline.size() == 0) {
            throw putils::FormatError(__FUNCTION__,name," blank line");
        }
        if (sline[0] == '#') {
            continue;
        }
        size_t ntokens = putils::tokenize_string(sline,delims,tokens);
        switch (ntokens)
        {
        case 0:
            throw putils::FormatError(__FUNCTION__,name,"not enough tokens on line");
        case 1:
            setValue(tokens[0],std::string("true"));
            break;
        case 2:
            setValue(tokens[0],tokens[1]);
            break;
        default:
            throw putils::FormatError(__FUNCTION__,name,"too many tokens on line");
        }
    }
    in.close();
}

void ProgramOptions::parseEnv( const char * prefix)
{
    std::string env_str("");
    std::string prefix_("");
    if ( prefix ) {
        prefix_ = std::string(prefix);
        prefix_ += "_";
    }
    std::cerr << "prefix " << prefix_ << "\n";
    for ( Option& opt : options)
    {
        env_str.clear();
        std::string name_ = opt.getName();
        std::string name = putils::string_toupper(name_);
        std::cerr << "name = " << name << "\n";
        env_str = prefix_;
        env_str += name;   
        std::cerr << "env_str = " << env_str << "\n";
        char * ret = getenv(env_str.c_str());
        if (ret) {
            std::cerr << "return = " << ret << "\n";
            opt.setValue(ret);
//            free(ret);
        }
    }
}

void ProgramOptions::parseCommandLine(int argc,char **argv)
{
    // look for null pointers
    if (argc) {
        if (argv == 0x0) {
            std::cerr << "null pointer argv\n";
            throw FormatError(__FUNCTION__,"null pointer");
        }
        for (int i=0;i<argc;++i) {
            if ( argv[i] == nullptr) {
                std::cerr << "null pointer argv[" << i << "]\n";
                throw FormatError(__FUNCTION__,"null pointer");
            }
            std::size_t sz = strlen(argv[i]);
            if (sz < 2) {
                std::cerr << "argv[" << i << "] is too short\n";
                std::cerr << argv[i] << "\n";
                throw FormatError(__FUNCTION__,"argument is too short");
            }
            if ( argv[i][1] == '-' && sz < 3) {
                std::cerr << "argv[" << i << "] is too short\n";
                std::cerr << argv[i] << "\n";
                throw FormatError(__FUNCTION__,"argument is too short");                
            }
        }
    }else{
        std::cerr << "argc = 0\n";
        throw FormatError(__FUNCTION__,"no arguments at all");
    }
    prog_name = std::string(argv[0]);
        for (auto k=1;k<argc;++k) {
            std::string str(argv[k]);
            if (str[0]=='-') {
                std::size_t p = (str[1] == '-') ? 2:1;
                std::string keyval = str.substr(p);
                if (keyval.compare("help")==0) printHelp();
                p = keyval.find("=");
                if (p!=std::string::npos) {
                    std::string key = keyval.substr(0,p);
                    std::string val = keyval.substr(p+1,std::string::npos);
                    setValue(key,val);
                } else {
                    if (k==(argc-1)) {
                        setValue(keyval,"true");
                    } else {
                        auto j = k + 1;
                        if (argv[j][0]=='-') {
                            setValue(keyval,"true");
                        } else {
                            ++k;
                            std::string val = std::string(argv[k]);
                            setValue(keyval,std::string(argv[k]));
                        }
                    }
                }
            }else{
                if (allow_extra_arguments) continue;
                std::cerr << "non option found on command line\n";
                throw FormatError(__FUNCTION__,"bad option");
            }
        }
}
} // end namespace putils
