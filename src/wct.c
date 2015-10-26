#include <string.h>
#include <math.h>
#include <assert.h>
#include "wct.h"
#include "wctparms.h"

int debug = 0;

/*Information about debug*/
int dbg_lvl()
{
    return debug;
}
void set_dbg_lvl( int dbglvl )
{
    debug = dbglvl;
}

/**
 * Reading of the jobfile
 */
static int permute_nodes ( int *invorder, int vcount, int *duration,
                           int *weight, int **durationlist, int **weightlist );
int read_problem( char *f, int *njobs, int **durationlist, int **weightlist )
{
    int i,val = 0;
    int nbjobs = 0,  prob = 0;
    int curduration, curweight, curjob = 0;
    int *duration = ( int * ) NULL;
    int *weight = ( int * ) NULL;
    double *ratio = ( double *) NULL;
    char buf[256], *p;
    int bufsize;
    const char *delim = " \n";
    char *data = ( char * ) NULL;
    char *buf2 = ( char * ) NULL;
    FILE *in = ( FILE * ) NULL;
    int *perm = ( int * ) NULL;
    int *iperm = ( int * )NULL;
    in = fopen( f, "r" );

    if ( !in ) {
        fprintf( stderr, "Unable to open file %s\n", f );
        val = 1;
        goto CLEAN;
    }

    if ( fgets( buf, 254, in ) != NULL ) {
        p = buf;

        if ( p[0] == 'p' ) {
            if ( prob ) {
                fprintf( stderr, "ERROR: in this file we have to p lines\n" );
                val = 1;
                goto CLEAN;
            }

            prob = 1;
            strtok( p, delim );
            data = strtok( NULL, delim );
            sscanf( data, "%d", &nbjobs );
            bufsize = 2 * nbjobs * ( 2 + ( int ) ceil( log( ( double )nbjobs + 10 ) ) );
            buf2 = ( char * ) CC_SAFE_MALLOC( bufsize, char );
            CCcheck_NULL_2( buf2, "Failed to allocate buf2" );
            weight = CC_SAFE_MALLOC( nbjobs, int );
            CCcheck_NULL_2( weight, "out of memory for weight" );
            duration = CC_SAFE_MALLOC(nbjobs, int);
            CCcheck_NULL_2(duration, "Failed to allocate memory");
            ratio = CC_SAFE_MALLOC(nbjobs, double);
            CCcheck_NULL_2(ratio, "Failed to allocate memory");
        } else {
            fprintf( stderr, "File has to give first the number vertices and edges.\n" );
            val = 1;
            goto CLEAN;
        }
    } else {
      val = 1;
      goto CLEAN;
    }

    while ( fgets( buf2, bufsize, in ) != ( char * )NULL ) {
        p = buf2;

        if ( p[0] == 'p' ) {
            if ( prob ) {
                fprintf( stderr, "ERROR: in this file we have to p lines\n" );
                val = 1;
                goto CLEAN;
            }
        } else
            if ( p[0] == 'n' ) {
                if ( !prob ) {
                    fprintf( stderr, "ERROR n before p in file\n" );
                    val = 1;
                    goto CLEAN;
                }

                strtok( p, delim );
                data = strtok( NULL, delim );
                sscanf( data, "%d", &curweight);
                data = strtok( NULL, delim );
                sscanf( data, "%d", &curduration);
                weight[curjob] = curweight;
                duration[curjob] = curduration;
                data = strtok( NULL, delim );
                curjob++;

            }
    }

    perm = CC_SAFE_MALLOC( nbjobs, int );
    CCcheck_NULL_2( perm, "Failed to allocate memory" );

    for (  i = 0; i < nbjobs; i++ ) {
        perm[i] = i;
        ratio[i] = (double) duration[i]/(double) weight[i];
    }
    


    CCutil_double_perm_quicksort(perm, ratio, nbjobs);
    iperm = CC_SAFE_MALLOC( nbjobs, int );
    CCcheck_NULL_2( iperm, "Failed to allocate memory" );

    for (  i = 0; i < nbjobs; ++i ) {
        iperm[perm[i]] = i;

    }

    permute_nodes( iperm, nbjobs, duration, weight, durationlist, weightlist );
    *njobs = nbjobs;

CLEAN:

    if ( val ) {
        CC_IFFREE( *durationlist, int );
        CC_IFFREE( *weightlist, int );
    }

    CC_IFFREE( weight, int );
    CC_IFFREE(duration, int);
    CC_IFFREE( buf2, char );
    CC_IFFREE( perm, int );
    CC_IFFREE( iperm, int );
    CC_IFFREE(ratio, double);

    if ( in ) {
        fclose( in );
    }

    return val;
}

static int permute_nodes ( int *invorder, int njobs,  int *duration,
                           int *weights, int **durationlist, int **weightlist )
{
    int i, val = 0;
    int *idurationlist = ( int * ) NULL, *iweightlist = ( int * ) NULL;

    iweightlist = CC_SAFE_MALLOC ( njobs, int );
    CCcheck_NULL_2 ( iweightlist, "out of memory for iweights" );
    idurationlist = CC_SAFE_MALLOC ( njobs, int );
    CCcheck_NULL_2 ( iweightlist, "out of memory for iweights" );


    for ( i = 0; i < njobs; i++ ) {
        iweightlist[invorder[i]] = weights[i];
        idurationlist[invorder[i]] = duration[i];
    }

    *durationlist =  idurationlist;
    *weightlist = iweightlist;
CLEAN:

    if ( val ) {
        CC_IFFREE ( iweightlist, int );
        CC_IFFREE ( idurationlist, int );
    }

    return val;
}

/*Functions for initialization of the problem and freeing the problem*/
void wctproblem_init( wctproblem *problem )
{
    /*B&B info*/
    problem->nwctdata = 0;
    problem-> global_upper_bound = INT_MAX;
    problem->first_lower_bound = 0;
    problem->global_lower_bound = 0;
    problem->first_upper_bound = INT_MAX;
    problem->status = no_sol;
    problem->rel_error = 1.0;
    problem->nbestschedule = 0;
    problem->bestschedule = ( Scheduleset *)NULL;
    problem->maxdepth = 0;
    problem->nbinitsets = 0;
    problem->gallocated = 0;
    problem->initsets = ( Scheduleset *)NULL;
    /*data of the problem*/
    wctdata_init( &( problem->root_pd ) );
    /*parms of the problem*/
    wctparms_init( &( problem->parms ) );
    
    /*heap initialization*/
    problem->br_heap = ( pmcheap *) NULL;
    pmcheap_init( &( problem->br_heap ), 1000 );
    problem->found = 0;
    
    /*CPU timer initialisation*/
    CCutil_init_timer( &( problem->tot_branching_cputime ),
                       "tot_branching_cputime" );
    CCutil_init_timer( &( problem->tot_cputime ), "tot_cputime" );
    CCutil_init_timer( &( problem->tot_lb_cpu_time ), "tot_lb_cpu_time" );
    CCutil_init_timer( &( problem->tot_lb ), "tot_lb" );
    CCutil_init_timer( &( problem->tot_scatter_search ), "tot_scatter_search" );
    CCutil_init_timer( &( problem->tot_pricing ), "tot_pricing" );
    
    /* Scatter sear init*/
    SS_init( &problem->scatter_search, 15, 10, 0.1 );
}

void wctproblem_free( wctproblem *problem )
{
    /*free the parameters*/
    wctparms_free( &( problem->parms ) );
    wctdata_free( &( problem->root_pd ) );
    /*free the heap*/
    pmcheap_free( problem->br_heap );
    problem->br_heap = ( pmcheap *) NULL;
    Schedulesets_free( &( problem->initsets ), &( problem->gallocated ) );
    Schedulesets_free( &( problem->bestschedule ), &( problem->nbestschedule ) );
    SS_free( &problem->scatter_search );
}

/*Functions for initialization and free the data*/
void wctdata_init( wctdata *pd )
{
    /*Initialization B&B data*/
    pd->id = -1;
    pd->depth = 0;
    pd->status = initialized;
    sprintf( pd->pname, "temporary" );
    /*Initialization graph data*/
    pd->njobs = 0;
    pd->orig_node_ids = ( int *) NULL;
    
    pd->duration = ( int *)NULL;
    pd->weights = ( int *)NULL;
    pd->upper_bound = INT_MAX;
    pd->lower_bound = 1;
    pd->dbl_safe_lower_bound = 0.0;
    pd->dbl_est_lower_bound = 0.0;
    pd->lower_scaled_bound = 1;
    pd->kpc_pi_scalef = 1;
    pd->LP_lower_bound = 0.0;
    /*Initialization  of the LP*/
    pd->LP = ( wctlp *)NULL;
    pd->x = ( double *)NULL;
    pd->coef = ( double *) NULL;
    pd->pi = ( double *) NULL;
    pd->kpc_pi = ( int *)NULL;
    /*Initialization pricing_problem*/

    /*Initialization of Scheduleset*/
    pd->ccount = 0;
    pd->cclasses = ( Scheduleset *)NULL;
    pd->dzcount = 0;
    pd->gallocated = 0;
    pd->newsets = ( Scheduleset *)NULL;
    pd->nnewsets = 0;
    pd->bestcolors = ( Scheduleset *) NULL;
    pd->nbestcolors = 0;
    pd->debugcolors = ( Scheduleset *) NULL;
    pd->ndebugcolors = 0;
    pd->opt_track = 0;
    /*Initialization max and retirement age*/
    pd->maxiterations = 1000000;
    pd->retirementage = 1000000;
    /*initialization of branches*/
    pd->v1 = pd->v2 = -1;
    pd->parent = ( wctdata *) NULL;
    pd->same_children = ( wctdata *) NULL;
    pd->nsame = 0;
    pd->diff_children = ( wctdata *) NULL;
    pd->ndiff = 0;
    heur_init( pd );
}

void lpwctdata_free( wctdata *pd )
{
    if ( pd->LP )
    {
        wctlp_free( &( pd->LP ) );
    }

    if ( pd->coef )
    {
        free( pd->coef );
        pd->coef = ( double *)NULL;
    }

    if ( pd->pi )
    {
        free( pd->pi );
        pd->pi = ( double *) NULL;
    }

    if ( pd->x )
    {
        free( pd->x );
        pd->x = ( double *)NULL;
    }


    if ( pd->kpc_pi )
    {
        free( pd->kpc_pi );
        pd->kpc_pi = ( int *)NULL;
    }

    heur_free( pd );
    Schedulesets_free( &( pd->newsets ), &( pd->nnewsets ) );
    Schedulesets_free( &( pd->cclasses ), &( pd->gallocated ) );
    pd->ccount = 0;
}

void children_data_free( wctdata *pd )
{
    int i;

    for ( i = 0; i < pd->nsame; ++i )
    {
        wctdata_free( &( pd->same_children[i] ) );
    }

    CC_IFFREE( pd->same_children, wctdata );

    for ( i = 0; i < pd->ndiff; ++i )
    {
        wctdata_free( &( pd->diff_children[i] ) );
    }

    CC_IFFREE( pd->diff_children, wctdata );
    pd->nsame = pd->ndiff = 0;
}

void temporary_data_free( wctdata *pd )
{
    children_data_free( pd );
    lpwctdata_free( pd );
}

void wctdata_free( wctdata *pd )
{
    temporary_data_free( pd );
    CC_IFFREE( pd->duration, int );
    CC_IFFREE( pd->orig_node_ids, int );
    CC_IFFREE( pd->weights, int );
    Schedulesets_free( &( pd->bestcolors ), &( pd->nbestcolors ) );
}

/** Help function for column generation */
void make_pi_feasible( wctdata *pd )
{
    int c;

    for ( c = 0; c < pd->ccount; ++c )
    {
        int i;
        double colsum = .0;
        double newcolsum = .0;

        for ( i = 0; i < pd->cclasses[c].count; ++i )
        {
            if ( signbit( pd->pi[pd->cclasses[c].members[i]] ) )
            {
                pd->pi[pd->cclasses[c].members[i]] = 0.0;
            }

            colsum += pd->pi[pd->cclasses[c].members[i]];
            colsum = nextafter( colsum, DBL_MAX );
        }

        if ( colsum > 1.0 )
        {
            for ( i = 0; i < pd->cclasses[c].count; ++i )
            {
                pd->pi[pd->cclasses[c].members[i]] /= colsum;
                newcolsum += pd->pi[pd->cclasses[c].members[i]];
            }

            if ( dbg_lvl() > 1 )
            {
                printf( "Decreased column sum of %5d from  %30.20f to  %30.20f\n", c, colsum,
                        newcolsum );
            }
        }
    }
}

static void reset_ages( Scheduleset *cclasses, int cccount )
{
    int i;

    for ( i = 0; i < cccount; i++ )
    {
        cclasses[i].age = 0;
    }
}

int add_newsets( wctdata *pd )
{
    int val = 0;
    Scheduleset *tmpsets = ( Scheduleset *) NULL;
    int i;

    if ( pd->nnewsets == 0 )
    {
        return val;
    }

    reset_ages( pd->newsets, pd->nnewsets );

    if ( pd->ccount + pd->nnewsets > pd->gallocated )
    {
        pd->gallocated *= 3;
        tmpsets = CC_SAFE_MALLOC( pd->gallocated, Scheduleset );
        CCcheck_NULL_2( tmpsets, "Failed to allocate memory to tmpsets" );
        memcpy( tmpsets, pd->cclasses, pd->ccount * sizeof( Scheduleset ) );
        free( pd->cclasses );
        pd->cclasses = tmpsets;
        tmpsets = NULL;
    }

    memcpy( pd->cclasses + pd->ccount, pd->newsets,
            pd->nnewsets * sizeof( Scheduleset ) );
    pd->ccount += pd->nnewsets;

    for ( i = pd->ccount; i < pd->gallocated; i++ )
    {
        Scheduleset_init( pd->cclasses + i );
    }

CLEAN:

    if ( val )
    {
        CC_IFFREE( pd->cclasses, Scheduleset );
    }

    CC_IFFREE( pd->newsets, Scheduleset );
    pd->nnewsets = 0;
    return val;
}

int double2int( int *kpc_pi, int *scalef, const double *pi, int vcount )
{
    int    i;
    double max_dbl_nweight = -DBL_MAX;
    double max_prec_dbl = exp2( DBL_MANT_DIG - 1 );
    static const double max_mwiswt   = ( double ) INT_MAX;
    double dbl_scalef = CC_MIN( max_prec_dbl, max_mwiswt );
    dbl_scalef /= ( double ) vcount;

    for ( i = 0; i < vcount; ++i )
    {
        max_dbl_nweight =
            CC_MAX( max_dbl_nweight, pi[i] );
    }

    dbl_scalef /= CC_MAX( 1.0, max_dbl_nweight );
    dbl_scalef  = floor( dbl_scalef );
    *scalef  = ( int ) dbl_scalef;

    for ( i = 0; i < vcount; ++i )
    {
        double weight = pi[i] * dbl_scalef;
        assert( weight < ( double ) INT_MAX );
        kpc_pi[i] = ( int ) weight;
    }

    return 0;
}

