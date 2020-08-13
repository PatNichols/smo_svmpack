#ifndef SVM_UTILS_H
#define SVM_UTILS_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_LINE_SIZE 128
#define MAX_TOKENS 128

inline void * Malloc(size_t n)
{
    void * ptr = malloc(n);
    if (ptr) return ptr;
    fprintf(stderr,"could not allocate %llu bytes\n",n);
    exit(EXIT_FAILURE);
}

inline void * Calloc(size_t n)
{
    void * ptr = malloc(n);
    if (ptr) {
        memset(ptr,0x0,n);
        return ptr;
    }
    fprintf(stderr,"could not allocate %llu bytes\n",n);
    exit(EXIT_FAILURE);
}

inline FILE *Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if (fp) return fp;
    fprintf(stderr,"could not open the file %s in mode %s\n",name,mode);
    exit(EXIT_FAILURE);
}


#define MALLOC_PTR(name) ( name *)Malloc(sizeof( name ) )

inline char ** tokens_init() {
    int i;
    char ** tok = (char**)Malloc(MAX_TOKENS*sizeof(char*));
    for (i=0; i<MAX_TOKENS; ++i) tok[i] = (char*)Malloc(MAX_LINE_SIZE);
    return tok;
}

inline void tokens_free(char **tokens)
{
    int i;
    for (i=MAX_TOKENS; i>0;) {
        --i;
        free(tokens[i]);
    }
    free(tokens);
}

inline int find_first_of(const char *str,const char *delims,int p)
{
    int i,j;
    char ch;
    int slen = strlen(str);
    int dlen = strlen(delims);

    for (i=p; i<slen; ++i) {
        ch = str[i];
        for (j=0; j<dlen; ++j) {
            if (ch==delims[j]) {
                return i;
            }
        }
    }
    return -1;
}

inline int find_first_not_of(const char *str,const char *delims,int p)
{
    int i,j;
    int is_inc;
    char ch;
    int slen = strlen(str);
    int dlen = strlen(delims);

    for (i=p; i<slen; ++i) {
        ch = str[i];
        is_inc = 0;
        for (j=0; j<dlen; ++j) {
            if (ch==delims[j]) {
                is_inc = 1;
                break;
            }
        }
        if (is_inc==0) return i;
    }
    return -1;
}

inline int explode_string(const char *str,const char *delims,char **tokens)
{
    int st,fn;
    int cnt= 0;
    size_t slen = strlen(str);
    st = find_first_not_of(str,delims,0);
    if (st==-1) return cnt;
    fn = find_first_of(str,delims,st);
    while (st!=-1 &&  cnt!=MAX_TOKENS) {
        if (fn==-1) {
            memcpy(tokens[cnt],str+st,strlen(slen-st);
                   tokens[cnt][slen]=0x0;
                   return cnt+1;
        }
        memcpy(tokens[cnt],str+st,fn-st);
        tokens[cnt][fn-st]=0x0;
        st = find_first_not_of(str,delims,fn);
        fn = find_first_of(str,delims,st);
    }
    return cnt+1;
}

typedef struct popt_t {
    char * key;
    char * des;
    char * val;
    int st;
};

typedef struct popt_list_node {
    popt_t * node;
    popt_t * next;
    popt_t * prev;
};

typedef struct popt_linked_list {
    popt_list_node * beg;
    popt_list_node * end;
};

typedef struct program_options_t {
    popt_linked_list * list;
    char * prog;
};

popt_t * popt_init(char *keyword,char *description,char *value)
{
    popt_t * p = MALLOC_PTR(popt_t);

    if (!key) {
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

inline int popt_matches(popt_t * p, char *flag) {
    return (memcmp(p->key,flag)==0);
}

inline popt_has_value(popt_t *p) {
    return p->st != -1;
}

inline int popt_was_set(popt_t *p) {
    return p->st == 0;
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


program_options_t * program_options_init(int argc,char **argv)
{
    program_options_t * opts = MALLOC_PTR(program_options_t);
    opts-> list = MALLOC_PTR(popt_linked_list);
    popt_linkd_list * list = (popt_linked_list*)Malloc(sizeof(popt_linked_list));
    list->beg = 0x0;
    list->end = 0x0;
    opts->prog = strdup(argv[0]);
}

void program_options_insert(program_options_t *opts,char *flag,char *des,char *val)
{
    popt_list_node * curr = opts->list->end;
    popt_list_node * p_node;
    if (!curr) {
        list->beg = (popt_list_node*)Malloc(sizeof(popt_list_node*));
        list->beg->prev = 0x0;
        list->beg->node = popt_init(flag,des,val);
        list->beg->next = 0x0;
        list->end = list->beg;
        return;
    }
    p_node = curr;
    curr->next = MALLOC_PTR(popt_list_node);
    curr = curr->next;
    curr->prev = p_node;
    curr->node = popt_init(flag,des,val);
    curr->next=0x0;
    list->end = curr;
    return;
}

void program_options_free(program_options *opts) {
    popt_list_node * curr = opts->list->end;
    while (curr) {
        p_node = curr;
        curr = curr->prev;
        popt_free(p_node->node);
        free(p_node);
    }
    free(opts->list);
    free(opts->prog);
    free(opts);
    return;
}

void program_options_write(program_options *opts,FILE *fp)
{
    popt_list_node * curr = opts->list->beg;
    fprintf(fp,"Usage is %s (options)\n",opts->prog);
    fprintf(fp,"Options are :\n");
    while (curr) {
        popt_write(curr->node,fp);
        curr = curr->next;
    }
}

popt_t * program_options_find(program_options_t *opt,char *flag)
{
    popt_list_node * curr = opts->list->beg;
    while (curr) {
        if (popt_matches(curr->node,flag)) return curr->node;
        curr = curr->next;
    }
    fprintf(stderr,"Could not find the option %s\n",flag);
    exit(EXIT_FAILURE);
}


popt_t * program_options_set_value(program_options_t *opt,char *flag,char *value)
{
    popt_list_node * curr = opts->list->beg;
    while (curr) {
        if (popt_matches(curr->node,flag)) {
            popt_set_value(curr->node,value);
            return curr->node;
        }
        curr = curr->next;
    }
    fprintf(stderr,"Could not find the option %s\n",flag);
    exit(EXIT_FAILURE);
}

char * program_options_get_value(program_options_t *opt,char *flag,char *value)
{
    popt_list_node * curr = opts->list->beg;
    while (curr) {
        if (popt_matches(curr->node,flag)) {
            return strdup(curr->node->val);
        }
        curr = curr->next;
    }
    fprintf(stderr,"Could not find the option %s\n",flag);
    exit(EXIT_FAILURE);
}

int program_options_has_value(program_options_t *opt,char *flag)
{
    popt_list_node * curr = opts->list->beg;
    while (curr) {
        if (popt_matches(curr->node,flag)) {
            return 1;
        }
        curr = curr->next;
    }
    return 0;
}

int program_options_was_set(program_options_t *opt,char *flag)
{
    popt_list_node * curr = opts->list->beg;
    while (curr) {
        if (popt_matches(curr->node,flag)) {
            return (curr->node->st==1);
        }
        curr = curr->next;
    }
    return 0;
}

inline void program_options_parse_command_line(program_options_t *opts,int argc,char **argv)
{
    int i,j;
    char key[256];
    char val[256];
    char * key_s;

    int eq_pos,pos,k,j,i;

    opts->prog = argv[0];
    for (i=0; i<argc; ++i) {
        key_s = argv[i];
        slen = strlen(key_s);
        if (slen > 1) {
            if (key_s[0]!='-') {
                parse_error("program options")
            }
            pos = (key_s[1]=='-') ? 2:1;
        }
        eq_pos = find(key_s+pos,'=');
        if (eq_pos==-1) {
            if (i!=(argc-1)) {
                if (argv[i+1][0]=='-') {
                    program_options_set_value(opts,key_s+pos,"true");
                } else {
                    program_options_set_value(opts,key_s+pos,argv[i+1]);
                }
            } else {
                program_options_set_value(opts,key_s+pos,"true");
            }
        } else {
            j = 0;
            for (k=pos; k<eq_pos; ++k,++j) key[j]=key_s[k];
            key[j] = 0x0;
            j = 0;
            for (k=eq_pos+1; k<slen; ++k,++j) val[k]=key_s[k];
            val[j]=0x0;
            program_options_set_value(opts,key,val);
        }
    }
}

inline void program_options_parse_config_file(program_options *opts,char *filename)
{
    const char *delims = " =:\n";
    char *buffer = (char*)Malloc(MAX_LINE_SIZE);
    char ** tokens;
    int nrd,ntokens;

    tokens_init(tokens);
    FILE *fp = Fopen(filename,"r");
    while (!feof(fp)) {
        nrd = getline(buffer,buffer_sz,fp);
        if (nrd==0) continue;
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
    return;
}

inline void program_options_parse_environment(program_options* opts,char *prefix)
{
    popt_list_node * curr = opts->list->beg;
    char flag[128];
    char *env_str;
    char *env_val;

    if (prefix==0x0) {
        while (curr) {
            env_str = strdup(curr->node->key);
            for (i=0; i<strlen(flag); ++i) {
                env_str[i] = toupper(env_str[i]);
            }
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(curr->node,env_val);
            }
            free(env_str);
            curr = curr->next;
        }
    } else {
        while (curr) {
            flag = strdup(curr->node->key);
            for (i=0; i<strlen(flag); ++i) {
                flag[i] = toupper(flag[i]);
            }
            env_str = strcat(prefix,"_");
            env_str = strcat(env_str,flag);
            env_val = getenv(env_str);
            if (env_val) {
                program_options_set_value(curr->node,env_val);
            }
            free(flag);
            curr = curr->next;
        }
    }
}

#endif