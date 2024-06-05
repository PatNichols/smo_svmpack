#pragma once
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <vector>
#include <cstdio>
#include "putils_c.h"
#include "putils_cxx.hpp"

namespace putils
{
struct NoSuchOption: public std::exception
{
    std::string msg;
    NoSuchOption() = delete;
    explicit NoSuchOption(const char * name):
        msg("no such option as ")
    {
        msg += name;
    }    
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

struct MissingValue: public std::exception
{
    std::string msg;
    MissingValue() = delete;
    explicit MissingValue(const char * name):
        msg("option ")
    {
        msg += name;
        msg += " does not have a value";
    }    
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

class ProgramOptions
{
    class Option
    {
        std::string m_name;
        std::string m_desc;
        std::string m_value;
        int m_status;
        int m_required;
    public:
        Option(const char *name, const char *description,const char *value,
                        int required = 0):
                m_name(name),
                m_desc( description ? description:std::string()),
                m_value( value ? value:std::string()),
                m_status( value ? 0:-1),
                m_required(required)
        {
                if ( m_value.size() ) m_status = 0;
        }

        Option(const std::string& name, const std::string& description,const std::string& value,
                        int required = 0):
                m_name(name),
                m_desc(description),
                m_value(value),
                m_status(-1),
                m_required(required)
        {
                if ( m_value.size() ) m_status = 0;
        }

        bool matches(const char *name) const noexcept
        { 
            return m_name.compare(name)==0;
        }

        bool matches(std::string& name) const noexcept
        { 
            return m_name.compare(name)==0;
        }

        constexpr bool required() const noexcept
        { 
            return m_required==1;
        }      

        void setValue(const char *value)
        {
//            std::cerr << "opt set value " << value << "\n";
            if (m_status != 1) {
                m_value = std::string(value);
//                std::cerr << " set value = " << m_value << "\n";
                m_status = 1;
            }
        }

        const std::string& getValue() const noexcept
        { 
            return m_value;
        }

        bool hasValue() const noexcept
        { 
            return m_value.size()!=0;
        }

        const std::string getName() const noexcept { return m_name;}

        std::ostream& write(std::ostream& os) const
        {
            os << m_name;
            if (m_desc.size()) {
                 os << " \"" << m_desc << "\"";
            }
            if (m_required) os << " required";
            os << "\n";
            os << "   ";
            if ( m_status != -1) {
                os << m_value;
                if ( m_status == 0) {
                    os << " default value";
                }else{
                    os << " set by user";
                }
            }else{
                os << " no value set";
            }
            os << "\n";
            return os;
        }

        friend std::ostream& operator << ( std::ostream& os, const Option& o)
        {
            return o.write(os);
        }
    };
    std::vector<Option> options;
    std::string prog_name;
    int allow_extra_arguments;
    int allow_extra_options;
public:
    ProgramOptions():
        options(),
        prog_name(),
        allow_extra_arguments(0),
        allow_extra_options(0)
    {}
    void setValue(const char *name,const char *value)
    {
        for ( Option& opt : options) {
            if ( opt.matches(name) ) 
            {
                opt.setValue(value);
                return;
            }
        }
        if ( allow_extra_options ) return;
        std::cerr << " no option " << name << "\n";
        throw NoSuchOption(name);
    }
    void setValue(const std::string& name,const std::string& value)
    {
        for ( Option& opt : options) {
            if ( opt.matches(name.c_str()) )
            {
                opt.setValue(value.c_str());
                return;
            }
        }
        if ( allow_extra_options ) return;
        std::cerr << " no option " << name << "\n";
        throw NoSuchOption(name.c_str());
    }

    bool hasRequiredValues() const noexcept
    {
        for (const Option& opt: options)
        {
            if ( opt.required() !=0 && !opt.hasValue())
            {
                return false;
            }
        }
        return true;
    }

    std::string getValue(const char *name)
    {
        for ( const Option& opt : options) {
            if ( opt.matches(name) ) {
                if ( opt.hasValue() ) 
                {
                    return opt.getValue();
                } 
                else 
                {
                    throw MissingValue(name);
                }
            }
        }
        throw NoSuchOption(name);
    }

    bool hasOption(const char *name) const
    {
        for ( const Option& o : options) {
            if (o.matches(name)) return true; 
        }
        return false;
    }

    bool hasValue(const char *name)
    {
        for ( const Option& o : options) {
            if (o.matches(name)) return o.hasValue(); 
        }
        throw NoSuchOption(name);
        return false;    
    }
    
    void addOption(const char *name, const char *description, const char *value,int required=0)
    {
        for (const Option& o: options) {
            if (o.matches(name)) {
                std::cerr << "warning : duplicate option!\n";
                std::cerr << o;
                std::cerr << "new option ";
                std::cerr << "name = " << name << "\n";
                if (description) {
                    std::cerr << std::string(description) << "\n";
                }
                if ( value ) std::cerr << value << "\n";
                return;
            }
        }
        options.push_back(Option(name,description,value,required));
    }
    void addOption(const std::string& name, const std::string& description, 
                    const std::string& value,int required=0)
    {
        for (const Option& o: options) {
            if (o.matches(name.c_str())) {
                std::cerr << "warning : duplicate option!\n";
                std::cerr << o;
                std::cerr << "new option ";
                std::cerr << "name = " << name << "\n";
                if (description.size()) {
                    std::cerr << std::string(description) << "\n";
                }
                if ( value.size() ) std::cerr << value << "\n";
                return;
            }
        }
        options.push_back(Option(name,description,value,required));
    }

    void addOptionsFromFile( const char *filename);    
    void parseCommandLine(int argc,char **argv);
    void parseFile( const char *file_name);
    void parseEnv( const char * prefix = nullptr);
    
    void setAllowExtraArguments(bool b=true) noexcept { allow_extra_arguments = (b) ? 1:0;}
    void setAllowExtraOptions(bool b=true) noexcept { allow_extra_options = (b) ? 1:0;}

    std::ostream& write(std::ostream& os) const
    {
        os << "usage is " << prog_name << " [OPTIONS]\n";
        for ( const Option& opt: options) os << opt;
        return os; 
    }

    void printHelp()
    {
        this->write(std::cerr);
        putils_error_function();
    }
    
    friend std::ostream& operator << ( std::ostream& os, const ProgramOptions& p)
    {
        return p.write(os);
    }
};
} // end namespace putils
