#ifndef PROGRAM_OPTIONS_H
#define PROGRAM_OPTIONS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"

typedef struct {
    char * key;
    char * des;
    char * val;
    int st;
} popt_t;

typedef struct {
    popt_t * v;
    size_t sz;
    size_t cap;
} popt_list_t;


typedef struct {
    popt_list_t * list;
    char * prog;
} program_options_t;


#define popt_matches(p,flag) strcmp((p)->key,(flag))==0

#define popt_has_value(p) (p)->st != -1

#define popt_was_set(p) (p)->st == 0

void popt_set_value(popt_t *p,char *new_value);
void popt_free(popt_t *p);
void popt_write(popt_t *p,FILE *fp);
popt_t * popt_init(char *keyword,char *description,char *value);
popt_t * popt_set(popt_t*p,char *,char *,char *);
popt_list_t * popt_list_allocate();
void popt_list_deallocate(popt_list_t *list);
void popt_list_grow(popt_list_t *list,size_t new_cap);
popt_t * popt_list_value(popt_list_t *list,size_t i);
void popt_list_insert(popt_list_t * list,char *k,char *d,char *v);

program_options_t * program_options_init(int argc,char **argv);
void program_options_insert(program_options_t *opts,char *flag,char *des,char *val);
void program_options_free(program_options_t *opts);
void program_options_write(program_options_t *opts,FILE *fp);
popt_t * program_options_find(program_options_t *opts,char *flag);
popt_t * program_options_set_value(program_options_t *opts,char *flag,char *value);
char * program_options_get_value(program_options_t *opts,char *flag);
int program_options_has_value(program_options_t *opts,char *flag);
int program_options_was_set(program_options_t *opts,char *flag);
void program_options_parse_command_line(program_options_t *opts,int argc,char **argv);
void program_options_parse_config_file(program_options_t* opts,char *filename);
void program_options_parse_environment(program_options_t* opts,char *prefix);
void program_options_insert_options(program_options_t*,char *);
void program_options_print_help(program_options_t *opts);

#endif