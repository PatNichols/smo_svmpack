#ifndef PUTILS_PROGRAM_OPTIONS_HPP
#define PUTILS_PROGRAM_OPTIONS_HPP
#include <vector>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <cctype>

namespace putils {

template < typename T >
inline T StringToType(const std::string& str) {
    try {
        std::istringstream is(str);
        T x;
        is >> x;
        return x;
    } catch (std::exception& e) {
        std::cerr << " StringToType cannot convert " << str << " to type requested\n";
        exit(EXIT_FAILURE);
    }
}

template <> inline int StringToType<int>(const std::string& str) {
    return std::stoi(str);
}
template <> inline unsigned int StringToType<unsigned int>(const std::string& str) {
    return static_cast<unsigned int>(std::stoi(str));
}
template <> inline long StringToType<long>(const std::string& str) {
    return std::stol(str);
}
template <> inline unsigned long StringToType<unsigned long>(const std::string& str) {
    return std::stoul(str);
}
template <> inline long long int StringToType<long long int>(const std::string& str) {
    return std::stoll(str);
}
template <> inline unsigned long long int StringToType<unsigned long long int>(const std::string& str)
{
    return std::stoull(str);
}

template <> inline float StringToType<float>(const std::string& str) {
    return std::stof(str);
}
template <> inline double StringToType<double>(const std::string& str) {
    return std::stod(str);
}
template <> inline long double StringToType<long double>(const std::string& str) {
    return std::stold(str);
}

template <> inline const char* StringToType< const char * >(const std::string& str) {
    return str.c_str();
}

template <> inline bool StringToType<bool>(const std::string& str) {
    if (str.size()==0) {
        std::cerr << "empty string passed to StringToType\n";
        exit(EXIT_FAILURE);
    }
    char ch = str[0];
    if (ch=='N' || ch=='n' || ch=='0' || ch=='F' || ch=='f') {
        return false;
    }
    return true;
}

template <> inline short StringToType<short>(const std::string& str) {
    int x = stoi(str);
    return (short)x;
}

template <> inline unsigned short StringToType<unsigned short>(const std::string& str) {
    int x = stoi(str);
    return (unsigned short)x;
}

template < typename T > inline std::string TypeToString( const T& x) {
    std::ostringstream out;
    out << x;
    return std::string(out.str());
}

inline void ToUpper(std::string& str) {
    for (auto j=0; j<str.size(); ++j) {
        if (isalpha(str[j])) str[j] = toupper(str[j]);
    }
}

inline void splitString(const std::string& str,const std::string& delims, std::vector<std::string>& tokens )
{
    tokens.clear();
    std::size_t st = str.find_first_not_of(delims,0);
    std::size_t fn = str.find_first_of(delims,st);
    while (st!=std::string::npos) {
        tokens.push_back(str.substr(st,fn-st));
        st = str.find_first_not_of(delims,fn);
        fn = str.find_first_of(delims,st);
    }
}

struct prog_opt {
    std::string key;
    std::string des;
    std::string val;
    int s;
    prog_opt(const std::string& key_str,const std::string& des_str):
        key(key_str),
        des(des_str),
        val(),s(-1) {}
    prog_opt(const std::string& key_str,const std::string& des_str, const std::string& val_str):
        key(key_str),
        des(des_str),
        val(val_str),s(0) {
        if (des.size()==0) des="no description";
        if (val.size()==0) s=-1;
    }

    bool matches(const std::string& flag) const noexcept {
        return (key.compare(flag)==0);
    }
    constexpr bool has_value() const noexcept {
        return s>=0;
    }
    constexpr bool was_set() const noexcept {
        return s==1;
    }
    void set_value(const std::string& val_str) noexcept {
        if (s==1) return;
        s=1;
        val = val_str;
    }
    inline const std::string& get_value() const noexcept {
        return val;
    }
    inline const std::string& GetDescription() const noexcept {
        return des;
    }
    inline const std::string& GetKeyword() const noexcept {
        return key;
    }
    inline std::ostream& WriteToStream(std::ostream& os) const {
        os << "-" << key << " : " << des << "\n     ";
        if (s==0) {
            os << val << " (default) ";
        }
        if (s==1) {
            os << val << " (set by user)";
        }
        os << "\n";
        return os;
    }
};


inline std::ostream& operator << (std::ostream& os,const prog_opt& opt) {
    return opt.WriteToStream(os);
}

struct NoSuchOptionError: public std::exception {
    const char *what() noexcept {
        return std::string("No Such Option was found\n").c_str();
    }
};


struct program_options {
    typedef std::vector< prog_opt >::iterator iterator;
    typedef std::vector< prog_opt >::const_iterator const_iterator;
    std::vector< prog_opt > opts;
    std::vector< prog_opt > unused_opts;
    std::vector< std::string > unused_words;
    std::string prog_name;
    bool allow_unused_args;
    bool allow_unused_opts;

    program_options():opts(),unused_opts(),unused_words(),prog_name("program_name"),
        allow_unused_args(false),allow_unused_opts(false) {}

    program_options(const std::string& opt_file):opts(),unused_opts(),unused_words(),prog_name("program_name"),
        allow_unused_args(false),allow_unused_opts(false) {
        add_options_file(opt_file);
    }

    void insert(const char *keyword,const char *descrip,const char *val) {
        opts.push_back(prog_opt(std::string(keyword),std::string(descrip),std::string(val)));
    }

    void AddOption(const std::string& keyword, const std::string& description = std::string("no description") ) {
        opts.push_back(prog_opt(keyword,description));
    }

    void AddOption(const std::string& keyword,const std::string& description,const std::string& value_str) {
        opts.push_back(prog_opt(keyword,description,value_str));
    }

    const std::string& get_value( const std::string& keyword) const noexcept {
        try {
            const_iterator iter = FindConstIterator(keyword);
            return iter->get_value();
        } catch (std::exception& e) {
            std::cerr << "program_options::get_value Could not find the option " << keyword << "\n";
            PrintHelp();
            exit(EXIT_FAILURE);
        }
    }

    bool has_value(const std::string& str) const noexcept {
        try {
            const_iterator iter = FindConstIterator(str);
            return iter->has_value();
        } catch (std::exception& e) {
            std::cerr << "Warning: program_options::has_value no option " << str << " exists\n";
            return false;
        }
    }

    void create_options_file() const {
        std::ofstream out("config.out");
        if (!out) {
            std::cerr << "Could not dump options!\n";
            exit(EXIT_FAILURE);
        }
        const_iterator iter = opts.cbegin();
        const_iterator iend = opts.end();
        while (iter!=iend) {
            out << iter->key << " " << iter->val <<"\n";
        }
        out.close();
        std::ofstream os("options.out");
        if (!os) {
            std::cerr << "Could not dump values for options!\n";
            exit(EXIT_FAILURE);
        }
        iter = opts.cbegin();
        while (iter!=iend) {
            os << iter->key << " \"" <<  iter->des << "\" " << iter->val <<"\n";
        }
        os.close();
        return;
    }


    void add_options_file(const std::string& OptionsFile)
    {
        std::ifstream in;
        try {
            in.open(OptionsFile.c_str(),std::ios::in);
        } catch (std::exception& e) {
            std::cerr << "Cannot read the file " << OptionsFile << " to insert options\n";
            exit(EXIT_FAILURE);
        }
        typedef std::size_t size_type;
        std::string sline;
        std::string delims(" =\n\0");
        std::vector< std::string > tokens;
        std::vector< std::string > new_tokens;

        while (in) {
            if (!getline(in,sline)) break;
            if (sline.size()==0) continue;
            if (sline[0]=='#') continue;
            splitString(sline,delims,tokens);
            size_type narg = tokens.size();
            size_type cnt = narg;
            // eliminate everything after a comment symbol
            for (size_type k=1; k<cnt; ++k) {
                if (tokens[k][0]=='#') {
                    narg = k;
                    break;
                }
            }
            if (narg!=cnt) {
                new_tokens.clear();
                for (size_type k=0; k<narg; ++k) {
                    new_tokens.push_back(tokens[k]);
                }
                tokens.swap(new_tokens);
                cnt = narg;
            }
            if (narg > 2 ) {
                char ch = tokens[1][0];
                if (ch=='\'' || ch=='\"') {
                    // doing collapse
                    ch = tokens[1][ tokens[1].size()-1 ];
                    if (ch=='\'' || ch=='\"') {
                        std::cerr << "warning description with quotes around a single word" << "\n";
                    } else {
                        cnt = 2;
                        new_tokens.clear();
                        new_tokens.push_back(tokens[0]);
                        std::string new_str(tokens[1]);
                        for (size_type k=2; k<narg; ++k) {
                            new_str += " ";
                            new_str += tokens[k];
                            ++cnt;
                            size_type m = tokens[k].size() - 1;
                            char ch = tokens[k][m];
                            if (ch=='\'' || ch=='\"') {
                                break;
                            }
                        }
                        new_tokens.push_back(new_str);
                        for (size_t k=cnt; k<narg; ++k) new_tokens.push_back(tokens[k]);
                        narg = new_tokens.size();
                        tokens.swap(new_tokens);
                    }
                }
            }
            switch (narg) {
            case 0:
                std::cerr << "warning empty line in parsing AddOptionFile \n" << sline << "\n";
                break;
            case 1:
                // assume we have only a keyword
                AddOption(tokens[0],"no description available");
                break;
            case 2:
                // we have a keyword + description
                AddOption(tokens[0],tokens[1]);
                break;
            case 3:
                // we have a keyword + description + default value
                AddOption(tokens[0],tokens[1],tokens[2]);
                break;
            default:
                std::cerr << "Error in input for AddOption File more than 3 tokens per line\n";
                exit(EXIT_FAILURE);
            }
            if (in.eof()) break;
        }
        in.close();
    }

    void parse_config_file(const std::string& ConfigFile)
    {
        std::ifstream in;
        try {
            in.open(ConfigFile.c_str(),std::ios::in);
        } catch (std::exception& e) {
            std::cerr << "ParseConfigFile:: Cannot read the file " << ConfigFile << " to insert options\n";
            exit(EXIT_FAILURE);
        }
        std::string sline;
        std::string delims(" =\n\0");
        std::vector<std::string> tokens;
        while (in) {
            if (!getline(in,sline)) break;
            std::size_t pos = sline.find_first_of("#");
            if (pos!=std::string::npos) {
                std::string sline2 = sline.substr(0,pos);
                sline.swap(sline2);
            }
            if (sline.size()==0) continue;
            splitString(sline,delims,tokens);
            std::size_t narg = tokens.size();
            switch (narg) {
            case 0:
                std::cerr << "ParseConfigFile : error reading line zero size\n" << sline << "\n";
                exit(EXIT_FAILURE);
            case 1:
                // assume we have a keyword and an assumed true value
                set_value(tokens[0],std::string("true"));
                break;
            case 2:
                // we have a keyword + value
                set_value(tokens[0],tokens[1]);
                break;
            default:
                std::cerr << "Error in input for ParseConfig File more than 2 tokens per line\n";
                std::cerr << sline << "\n";
                exit(EXIT_FAILURE);
            }
            if (in.eof()) break;
        }
        in.close();
    }

    void parse_environment(const std::string& EnvPrefix=std::string(""))
    {
        if (EnvPrefix.size()==0) {
            std::size_t nopts = opts.size();
            for (auto k=0; k<nopts; ++k) {
                std::string key = opts[k].key;
                putils::ToUpper(key);
                char * env_value = getenv(key.c_str());
                if (env_value) {
                    opts[k].set_value(std::string(env_value));
                }
            }
        } else {
            std::size_t nopts = opts.size();
            for (auto k=0; k<nopts; ++k) {
                std::string key = opts[k].key;
                putils::ToUpper(key);
                key = EnvPrefix + "_" + key;
                char * env_value = getenv(key.c_str());
                if (env_value) {
                    opts[k].set_value(std::string(env_value));
                }
            }
        }
    }

    void parse_command_line(int argc,char **argv)
    {
        prog_name = std::string(argv[0]);
        for (auto k=1; k<argc; ++k) {
            std::string str(argv[k]);
            if (str[0]=='-') {
                std::size_t p = (str[1] == '-') ? 2:1;
                std::string keyval = str.substr(p);
                if (keyval.compare("help")==0) PrintHelp();
                p = keyval.find("=");
                if (p!=std::string::npos) {
                    std::string key = keyval.substr(0,p);
                    std::string val = keyval.substr(p+1,std::string::npos);
                    set_value(key,val);
                } else {
                    if (k==(argc-1)) {
                        set_value(keyval,"true");
                    } else {
                        auto j = k + 1;
                        if (argv[j][0]=='-') {
                            set_value(keyval,"true");
                        } else {
                            ++k;
                            set_value(keyval,std::string(argv[k]));
                        }
                    }
                }
            } else {
                if (allow_unused_args) {
                    unused_words.push_back(str);
                    continue;
                }
                std::cerr << "Parse Error expected an -option_name but found a word not preceeded by a -\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    iterator FindIterator(const std::string& name) {
        iterator iter = opts.begin();
        while (iter!=opts.end()) {
            if (iter->matches(name)) return iter;
            ++iter;
        }
        std::string err_msg("no such option as ");
        err_msg += name;
        throw std::invalid_argument(err_msg);
    }

    const_iterator FindConstIterator(const std::string& name) const {
        const_iterator iter = opts.cbegin();
        while (iter!=opts.cend()) {
            if (iter->matches(name)) return iter;
            ++iter;
        }
        std::string err_msg("no such option as ");
        err_msg += name;
        throw std::invalid_argument(err_msg);
    }

    inline void set_value(const std::string& opt_name,const std::string& opt_value) noexcept {
        try {
            iterator iter = FindIterator(opt_name);
            iter->set_value(opt_value);
        } catch (std::exception& e) {
            if (allow_unused_opts) return;
            std::cerr << "program_options::set_value  no such option as " << opt_name <<"\n";
            PrintHelp();
        }
    }

    std::ostream& WriteToStream(std::ostream& os) const {
        os <<"Usage is :\n" << prog_name << " [ options ] \n";
        os <<"Options are :\n";
        const_iterator iter = opts.cbegin();
        while (iter!=opts.cend()) {
            os << *iter << "\n";
            ++iter;
        }
        return os;
    }

    void PrintHelp() const noexcept {
        this->WriteToStream(std::cerr);
        exit(EXIT_FAILURE);
    }

};

inline std::ostream& operator << ( std::ostream& os, const program_options& p) {
    return p.WriteToStream(os);
}

} // end namespace
#endif