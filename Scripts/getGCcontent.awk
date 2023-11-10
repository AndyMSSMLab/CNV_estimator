##Script to generate GC content per fasta sequence. 
awk ' \
BEGIN { \
    FS=""; \
    cg=0; \
    t=0; \
} \
{ \
    if ($1 != ">") { \
        for (i = 1; i <= NF; i++) { \
            if ($i ~ /[ACTGactg]/) { \
                t++;
            } \
            if ($i ~ /[CGcg]/) { \
                cg++;
            } \
        } \
    } \
    else { \
        if (t > 0) { \
            printf("%s\t%d\t%d\t%f\t%.0f\n",h,cg,t,cg/t,100*(cg/t))
            cg = 0; \
            t = 0; \
        } \
        h = substr($0,2); \
    } \
} \
END { \
    printf("%s\t%d\t%d\t%f\t%.0f\n",h,cg,t,cg/t,100*(cg/t))
}' $1
