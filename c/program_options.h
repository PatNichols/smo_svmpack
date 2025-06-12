#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"

struct popt_t {
    char * key;
    char * des;
    char * val;
    int st;
};

typedef struct popt_t popt_t;

struct popt_list_t{
    popt_t * v;
    size_t sz;
    size_t cap;
};

typedef struct popt_list_t popt_list_t;

struct program_options_t {
    popt_list_t * list;
    char * prog;
};

typedef struct program_options_t program_options_t;

program_options_t * program_options_init();
void program_options_free(program_options_t *opts);

void program_options_insert(program_options_t *opts,const char *flag,const char *des,const char *val);
void program_options_write(const program_options_t *opts,FILE *fp);
popt_t * program_options_find(program_options_t *opts,const char *flag);
void program_options_set_value(program_options_t *opts,const char *flag,const char *value);
const char * program_options_get_value(program_options_t *opts,const char *flag);
int program_options_has_value(const program_options_t *opts,const char *flag);
int program_options_was_set(const program_options_t *opts,const char *flag);
void program_options_parse_command_line(program_options_t *opts,int argc,char **argv);
void program_options_parse_config_file(program_options_t* opts,const char *filename);
void program_options_parse_environment(program_options_t* opts,const char *prefix);
void program_options_insert_options(program_options_t*,const char *);
void program_options_print_help(const program_options_t *opts);

#endif