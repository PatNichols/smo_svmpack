#include "program_options.h"

popt_t * popt_init(char *keyword,char *description,char *value)
{
    popt_t * p = MALLOC_PTR(popt_t);
    if (!keyword) {
        fprintf(stderr,"Cannot initialize a option with no keyword\n");
        exit(EXIT_FAILURE);
    }
    p-> key = strdup(keyword);
    if (description) {
        p-> des = strdup(description);
    } else {
        p-> des = 0x0;
    }
    if (value) {
        p-> val = strdup(value);
        p-> st = 0;
    } else {
        p->val = 0x0;
        p->st = -1;
    }
    return p;
}

inline popt_t * popt_set(popt_t *p,char *keyword,char *description,char *value)
{
    if (!keyword) {
        fprintf(stderr,"Cannot initialize a option with no keyword\n");
        exit(EXIT_FAILURE);
    }
    p-> key = strdup(keyword);
    if (description) {
        p-> des = strdup(description);
    } else {
        p-> des = 0x0;
    }
    if (value) {
        p-> val = strdup(value);
        p-> st = 0;
    } else {
        p->val = 0x0;
        p->st = -1;
    }
    return p;
}



inline void popt_set_value(popt_t *p,char *new_value) {
    if (p->st == 1) return;
    if (p->val) free(p->val);
    p->val = strdup(new_value);
    p->st = 1;
}

inline void popt_free(popt_t *p) {
    if (p->val) free(p->val);
    if (p->des) free(p->des);
    free(p->key);
    free(p);
}

inline void popt_write(popt_t *p,FILE *fp)
{
    fprintf(fp,"-%s :",p->key);
    if (p->des) {
        fprintf(fp,"%s\n",p->des);
    } else {
        fprintf(fp,"no description\n");
    }
    if (p->val) {
        fprintf(fp,"   value = %s",p->val);
        if (p->st == 0) {
            fprintf(fp,"(default)\n");
        } else {
            fprintf(fp,"(set by user)\n");
        }
    } else {
        fprintf(fp,"   no value\n");
    }
}

popt_list_t * popt_list_allocate() {
    popt_list_t * list = MALLOC_PTR(popt_list_t);
    list->cap = 60;
    list->v = (popt_t*)Malloc(sizeof(popt_t)*list->cap);
    list->sz = 0;
    return list;
}

void popt_list_deallocate(popt_list_t *list) {
    popt_t *p;
    size_t i;

    for (i=list->sz; i>0;) {
        i--;
        p = list->v + i;
        free(p->val);
        free(p->des);
        free(p->key);
    }
    free(list->v);
    free(list);
}

void popt_list_grow(popt_list_t *list,size_t new_cap) {
    popt_t * ptmp;
    popt_t * p;
    int i;
    if (new_cap == 0) new_cap = list->cap + 20;
    if (list->sz == list->cap) {
        ptmp = (popt_t*)Malloc(sizeof(popt_t)*new_cap);
        memcpy(ptmp,list->v,sizeof(popt_t)*list->sz);
        for (i=list->sz;i>0;) {
            i--;
            p = list->v + i;
            free(p->val);
            free(p->des);
            free(p->key);
        }
        free(list->v);
        list->v = ptmp;
    }
    list->cap = new_cap;
}

popt_t * popt_list_value(popt_list_t *list,size_t i) {
    return list->v+i;
}

void popt_list_insert(popt_list_t * list,char *k,char *d,char *v)
{
    popt_t * p;
    if (list->sz == list->cap) {
        popt_list_grow(list,20);
    }
    p = list->v + list->sz;
    popt_set(p,k,d,v);
    list->sz += 1;
}

inline void program_options_print_help(program_options_t *opts)
{
    program_options_write(opts,stderr);
    fprintf(stderr,"-help : print this screen\n");
    exit(EXIT_FAILURE);
}

inline void program_options_parse_command_line(program_options_t *opts,int argc,char **argv)
{
    int i,j;
    char *key;
    char *val;
    char * key_s;

    int eq_pos,pos,k,slen;

    opts->prog = strdup(argv[0]);
    fprintf(stderr,"program %s\n",argv[0]);
    for (i=1; i<argc; ++i) {
        key_s = argv[i];
        slen = strlen(key_s);
        if (slen > 1) {
            if (key_s[0]!='-') {
                parse_error("program options BAD OPTION no -");
            }
            pos = (key_s[1]=='-') ? 2:1;
            key_s = key_s  + pos;
            if (strcmp(key_s,"help")==0) program_options_print_help(opts);
            eq_pos = find_first_of(key_s,"=",0);
            if (eq_pos==-1) {
                if (i!=(argc-1)) {
                    if (argv[i+1][0]=='-') {
                        program_options_set_value(opts,key_s,"true");
                    } else {
                        program_options_set_value(opts,key_s,argv[i+1]);
                        ++i;
                    }
                } else {
                    program_options_set_value(opts,key_s,"true");
                }
            } else {
                key = strndup(key_s,eq_pos);
                val = strdup(key_s+eq_pos + 1);
                program_options_set_value(opts,key,val);
                free(val);
                free(key);
            }
        } else {
            parse_error("program options BAD ARGV");
        }
    }
}

inline void program_options_parse_config_file(program_options_t *opts,char *filename)
{
    const char *delims = " =:\n";
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char ** tokens;
    int nrd,ntokens;

    tokens = tokens_init();
    FILE *fp = Fopen(filename,"r");
    while (!feof(fp)) {
        nrd = getline(&buffer,&buffer_size,fp);
        if (nrd==0 || strlen(buffer)==0) break;
        if (buffer[0]=='#') continue;
        ntokens = explode_string(buffer,delims,tokens);
        if (ntokens==0) continue;
        if (ntokens!=2) {
            fprintf(stderr,"Error reading config file too many word on line\n");
            exit(EXIT_FAILURE);
        }
        program_options_set_value(opts,tokens[0],tokens[1]);
    }
    fclose(fp);
    tokens_free(tokens);
    free(buffer);
    fprintf(stderr,"parsed config file %s\n",filename);
    return;
}

inline void program_options_parse_environment(program_options_t* opts,char *prefix)
{
    char *flag;
    char *env_str;
    char *env_val;
    int i;
    size_t j;
    size_t nopt = opts->list->sz;
    popt_t * p = opts->list->v;

    if (prefix==0x0) {
        for (j=0; j<nopt; ++j) {
            env_str = strdup((p+j)->key);
            for (i=0; i<strlen(env_str); ++i) {
                env_str[i] = toupper(env_str[i]);
            }
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(opts,(p+j)->key,env_val);
                free(env_str);
            }
        }
    } else {
        for (j=0; j<nopt; ++j) {
            env_str = strcat(prefix,"_");
            env_str = strcat(env_str,(p+j)->key);
            for (i=0; i<strlen(env_str); ++i) {
                env_str[i] = toupper(env_str[i]);
            }
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(opts,(p+j)->key,env_val);
            }
        }
    }
}

inline void program_options_insert_options(program_options_t *popt,char *filename)
{
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char **tokens = tokens_init();
    int nrd,ntokens;
    FILE *fp = fopen(filename,"r");
    const char *delims = " =:\n";

    while (!feof(fp)) {
        nrd = getline(&buffer,&buffer_size,fp);
        if (nrd) {
            if (buffer[0]=='#') continue;
            ntokens = explode_string(buffer,delims,tokens);
            if (ntokens<=3) {
                switch(ntokens) {
                case 0:
                    break;
                case 1:
                    program_options_insert(popt,tokens[0],0x0,0x0);
                    break;
                case 2:
                    program_options_insert(popt,tokens[0],tokens[1],0x0);
                    break;
                case 3:
                    program_options_insert(popt,tokens[0],tokens[1],tokens[2]);
                    break;
                }
            } else {
                parse_error("program_options_insert_options");
            }
        }
    }
    fclose(fp);
}

inline program_options_t * program_options_init(int argc,char **argv)
{
    program_options_t * opts = MALLOC_PTR(program_options_t);
    opts->list = popt_list_allocate();
    opts->prog = strdup(argv[0]);
    return opts;
}

inline void program_options_insert(program_options_t *opts,char *flag,char *des,char *val)
{
    popt_list_insert(opts->list,flag,des,val);
    return;
}

inline void program_options_free(program_options_t *opts) {
    popt_list_deallocate(opts->list);
    free(opts->prog);
    free(opts);
    return;
}

inline void program_options_write(program_options_t *opts,FILE *fp)
{
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v = opts->list->v;
    fprintf(fp,"Usage is %s (options)\n",opts->prog);
    fprintf(fp,"Options are :\n");
    for (i=0; i<nopt; ++i) {
        fprintf(fp," name = %s\n",v[i].key);
        popt_write(&v[i],fp);
    }
}

inline popt_t * program_options_find(program_options_t *opts,char *flag)
{
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v= opts->list->v;
    for (i=0; i<nopt; ++i) {
        if (popt_matches(v+i,flag)) return v+i;
    }
    fprintf(stderr,"Could not find the option %s\n",flag);
    exit(EXIT_FAILURE);
}


inline popt_t * program_options_set_value(program_options_t *opts,char *flag,char *value)
{
    popt_t * p = program_options_find(opts,flag);
    popt_set_value(p,value);
    return p;
}

inline char * program_options_get_value(program_options_t *opts,char *flag)
{
    popt_t * p = program_options_find(opts,flag);
    return p->val;
}

inline int program_options_has_value(program_options_t *opts,char *flag)
{
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v= opts->list->v;
    for (i=0; i<nopt; ++i) {
        if (popt_matches(v+i,flag))
            return (v+i)->st >= 0;
    }
    fprintf(stderr,"warning no option %s was found\n",flag);
    return 0;
}

inline int program_options_was_set(program_options_t *opts,char *flag)
{
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v= opts->list->v;
    for (i=0; i<nopt; ++i) {
        if (popt_matches(v+i,flag))
            return (v+i)->st == 1;
    }
    fprintf(stderr,"warning no option %s\n",flag);
    return 0;
}
