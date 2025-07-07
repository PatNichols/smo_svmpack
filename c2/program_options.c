#include "program_options.h"

#define popt_matches(p,flag) strcmp((p)->key,(flag))==0

inline void popt_set_value(popt_t *p,const char *v) {
    if (p->st!=1) {
        p->st = 1;
        if (p->val) free(p->val);
        p->val = Strdup(v);
    }
}

inline void popt_free(popt_t *p) {
    if (p->val) free(p->val);
    if (p->des) free(p->des);
    if (p->key) free(p->key);
}

inline void popt_set(popt_t *p,const char *name,const char *descript,const char *value) {
    p->key = Strdup(name);
    p->des = (descript) ? Strdup(descript):0x0;
    if (value) {
        p->val = Strdup(value);
        p->st = 0;
    } else {
        p->val = 0x0;
        p->st = -1;
    }
}

inline void popt_clone(popt_t *new_p, const popt_t *p) {
    new_p->key = Strdup(p->key);
    new_p->des = (p->des) ? Strdup(p->des):0x0;
    new_p->val = (p->val) ? Strdup(p->val):0x0;
    new_p->st = p->st;
}

inline void popt_write(const popt_t *p,FILE *fp) {
    fprintf(fp,"- %s",p->key);
    if (p->des) fprintf(fp," %s",p->des);
    fprintf(fp,"\n");
    if (p->val) {
        fprintf(fp," value = %s",p->val);
        if (p->st == 0) fprintf(fp," default");
        else fprintf(fp," set by user");
    } else {
        fprintf(fp," no value set");
    }
    fprintf(fp,"\n");
}

inline void popt_list_deallocate(popt_list_t* list) {
    int i;
    popt_t * curr;
    for (i=(list)->cap; i;) {
        --i;
        curr = (list)->v + i;
        popt_free(curr);
    }
    free((list)->v);
}

inline void popt_list_grow(popt_list_t *list,size_t new_cap) {
    int i;
    popt_t * curr;
    popt_t * old_v = (list)->v;
    popt_t * new_v = (popt_t *)Malloc(sizeof(popt_t)*new_cap);
    for (i=0; i<(list)->sz; ++i) {
        curr = old_v+i;
        popt_clone(new_v+i,curr);
        popt_free(curr);
    }
    for (i=((list)->sz); i<new_cap; ++i) {
        (new_v+i)->key = 0x0;
        (new_v+i)->des = 0x0;
        (new_v+i)->val = 0x0;
    }
    (list)->cap = new_cap;
    (list)->v = new_v;
    free(old_v);
}


inline void popt_list_insert(popt_list_t * list,const char *key,const char *des,const char *val) {
    popt_t * opt;
    size_t new_cap;
    if ( list->sz >= list->cap) {
        new_cap = list->cap + list->cap;
        popt_list_grow(list,new_cap);
    }
    popt_set( (list->v + list->sz), key, des, val);
    list->sz += 1;
}

popt_list_t * popt_list_allocate() {
    size_t i;
    const size_t cap = 32UL;
    popt_t * curr;
    popt_list_t * list = MALLOC_PTR(popt_list_t);
    list->cap = cap;
    list->v = (popt_t*)Malloc(sizeof(popt_t)*cap);
    list->sz = 0;
    for (i=0; i<cap; ++i) {
        curr = list->v + i;
        curr->key = 0x0;
        curr->des = 0x0;
        curr->val = 0x0;
        curr->st = -1;
    }
    return list;
}

void program_options_print_help(const program_options_t *opts) {
    program_options_write(opts,stderr);
    fprintf(stderr,"-help : print this screen\n");
    exit(EXIT_FAILURE);
}

void program_options_parse_command_line(program_options_t *opts,int argc,char **argv) {
    int i,j;
    char *key;
    char *val;
    char * key_s;
    char * eq_ptr;
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
            eq_ptr = strstr(key_s,"=");
            if (eq_ptr == 0x0) {
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
                eq_pos = eq_ptr - key_s;
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

void program_options_parse_config_file(program_options_t *opts,const char *filename) {
    const char *delims = " \n";
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char ** tokens;
    int nrd,ntokens;
    char * ptr;
    tokens = tokens_init();
    FILE *fp = Fopen(filename,"r");
    while (!feof(fp)) {
        nrd = getline(&buffer,&buffer_size,fp);
        if (nrd==0 || strlen(buffer)==0) break;
        if (buffer[0]=='#') continue;
        ptr = strchr(buffer,'=');
        if (ptr) *ptr = ' ';
        ntokens = explode_string(buffer,delims,tokens);
        switch (ntokens) {
        case 0:
            fprintf(stderr,"format error : %s\n",buffer);
            exit(EXIT_FAILURE);
        case 1:
            program_options_set_value(opts,tokens[0],"true");
            break;
        case 2:
            program_options_set_value(opts,tokens[0],tokens[1]);
            break;
        default:
            fprintf(stderr,"format error : %s\n",buffer);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);
    tokens_free(tokens);
    free(buffer);
    fprintf(stderr,"parsed config file %s\n",filename);
    return;
}

void program_options_parse_environment(program_options_t* opts,const char *prefix) {
    char *flag;
    char *env_str;
    char *env_val;
    int i;
    size_t j;
    size_t nopt = opts->list->sz;
    popt_t * p = opts->list->v;
    char * prestr;
    
    if (prefix==0x0) {
        for (j=0; j<nopt; ++j) {
            env_str = Strdup((p+j)->key);
            for (i=0; i<strlen(env_str); ++i) {
                env_str[i] = toupper(env_str[i]);
            }
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(opts,(p+j)->key,env_val);
            }
            free(env_str);
        }
    } else {
        prestr = Strdup(prefix);
        env_str = (char *)Malloc(128);
        strcat(prestr,"_");
        for (j=0; j<nopt; ++j) {
            strcpy(env_str,prestr);
            strcat(env_str,p->key);
            for (i=0; i<strlen(env_str); ++i) {
                env_str[i] = toupper(env_str[i]);
            }
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(opts,(p+j)->key,env_val);
            }
        }
        free(env_str);
    }
}

void program_options_insert_options(program_options_t *popt,const char *filename) {
    size_t buffer_size = MAX_LINE_SIZE;
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char **tokens = tokens_init();
    char *desc = 0x0;
    char *beg = 0x0;
    char *end = 0x0;
    char *last = 0x0;
    int nrd,ntokens;
    FILE *fp = fopen(filename,"r");
    const char *delims = " \n";
    char quote_char = '\"';
    while (!feof(fp)) {
        nrd = getline(&buffer,&buffer_size,fp);
        if (nrd == 0 ) continue;
        if (buffer[0] == '#') continue;
        if (nrd > 0) {
            // look for desription string and then extract it from buffer
            beg = strchr(buffer,quote_char);
            if (beg) {
                end = strrchr(buffer,quote_char);
                if (beg == end) {
                    fclose(fp);
                    fprintf(stderr,"unterminated quoted string");
                    exit(EXIT_FAILURE);
                }
                last = Strdup(end+1);
                desc = strndup(beg+1,(end-beg-1));
                strcat(beg," ");
                strcat(beg+1,last);
                free(last);
            }
            ntokens = explode_string(buffer,delims,tokens);
            if (ntokens<=2) {
                switch(ntokens) {
                case 0:
                    break;
                case 1:
                    program_options_insert(popt,tokens[0],desc,0x0);
                    break;
                case 2:
                    program_options_insert(popt,tokens[0],desc,tokens[1]);
                default:
                    break;
                }
            } else {
                parse_error("program_options_insert_options");
            }
        } else {
            if ( nrd == -1) {
                if (feof(fp)) break;
                fprintf(stderr,"read error in %s",__FUNCTION__);
                exit(-1);
            }
        }
        if (desc) {
            free(desc);
            desc = 0x0;
        }
    }
    if (desc) free(desc);
    fclose(fp);
}

program_options_t * program_options_init() {
    program_options_t * opts = MALLOC_PTR(program_options_t);
    opts->list = popt_list_allocate();
    opts->prog = 0x0;
    return opts;
}

void program_options_insert(program_options_t *opts,const char *flag,
    const char *des,const char *val) {
    popt_list_insert(opts->list,flag,des,val);
    return;
}

void program_options_free(program_options_t *opts) {
    popt_list_deallocate(opts->list);
    if (opts->prog) free(opts->prog);
    free(opts);
    return;
}

void program_options_write(const program_options_t *opts,FILE *fp) {
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

popt_t * program_options_find(program_options_t *opts,const char *flag) {
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v= opts->list->v;
    for (i=0; i<nopt; ++i) {
        if (popt_matches(v+i,flag)) return v+i;
    }
    fprintf(stderr,"Could not find the option %s\n",flag);
    exit(EXIT_FAILURE);
}


void program_options_set_value(program_options_t *opts,
    const char *flag,const char *value) {
    popt_t * p = program_options_find(opts,flag);
    popt_set_value(p,value);
}

const char * program_options_get_value(program_options_t *opts,const char *flag) {
    const popt_t * p = program_options_find(opts,flag);
    return p->val;
}

int program_options_has_value(const program_options_t *opts,const char *flag) {
    size_t i;
    size_t nopt = opts->list->sz;
    popt_t * v= opts->list->v;
    for (i=0; i<nopt; ++i) {
        if (popt_matches(v+i,flag)) {
            return (v+i)->st >= 0;
        }
    }
    fprintf(stderr,"warning no option %s was found\n",flag);
    return -1;
}
