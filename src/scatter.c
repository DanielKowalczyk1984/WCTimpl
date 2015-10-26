#include <math.h>
#include <assert.h>
#include "wct.h"

static int select_parrent( GRand *Rand, int counter, int *sum, int totsum );
static int minimum_partition( solution *sol, int lowerbound );

int sort_vertices( const void *a, const void *b, void *data );

static int nodepair_ref_key( int v1, int v2 );
static void inodepair_ref_key( int *v1, int *v2, int index );

void print_totalweightcomptime( void *data, void *user_data );
int order_totalweightcomptime( const void *a, const void *b, void *data );
int order_totalweightcomptime_list( const void *a, const void *b );


/** compute row-index v1 and column-index v2 from array-index.*/
static void inodepair_ref_key( int *v1, int *v2, int index )
{
    *v2 = ( int ) floor ( sqrt( 2 * ( ( double )index ) + 0.25 ) - 0.5 );
    *v1 = index - ( *v2 * ( *v2 + 1 ) / 2 );
}

static int nodepair_ref_key( int v1, int v2 )
{
    /* We store only the elements of the upper right triangle within the
     njobs x njobs matrix. */
    assert( v1 <= v2 );
    return v2 * ( v2 + 1 ) / 2 + v1;
}

static int minimum_partition( solution *sol, int lowerbound )
{
    int val = 0;

    for ( int i = 1; i < sol->nmachines; ++i ) {
        if ( sol->part[i].completiontime - lowerbound <= 0 ) {
            if ( ( CC_OURABS( sol->part[i].completiontime - lowerbound ) <
                    CC_OURABS( sol->part[val].completiontime - lowerbound ) )
                    || ( ( CC_OURABS( sol->part[i].completiontime - lowerbound ) ==
                           CC_OURABS( sol->part[val].completiontime - lowerbound ) &&
                           g_queue_get_length( sol->part[i].list ) > g_queue_get_length(
                               sol->part[val].list ) ) ) ) {
                val = i;
            }
        }
    }
    return val;
}
/*
For each functions
 */

void distance_min_max( void *data, void *user_data )
{
    int k, l;
    min_max *temp = ( min_max * )user_data;
    solution *new_sol = temp->new_sol;
    int n = temp->n;
    solution *sol = ( solution * )data;
    int curval = 0;

    for ( l = 0; l < n - 1; l++ ) {
        for ( k = l + 1; k < n; k++ ) {
            if ( ( new_sol->vlist[l].part == new_sol->vlist[k].part
                    && sol->vlist[l].part != sol->vlist[k].part ) ||
                    ( new_sol->vlist[l].part != new_sol->vlist[k].part
                      && sol->vlist[l].part == sol->vlist[k].part ) ) {
                curval++;
            }

            if ( curval > temp->min ) {
                break;
            }
        }

        if ( curval > temp->min ) {
            break;
        }
    }

    if ( curval <= temp->min ) {
        temp->min = curval;
    }
}

void max_dist( void *data, void *user_data )
{
    solution *sol = ( solution * ) data;
    solution *max = ( solution * ) user_data;

    if ( ( sol->dist > max->dist ) || ( sol->dist == max->dist
                                        && sol->totalweightcomptime < max->totalweightcomptime ) ) {
        max = data;
    }
}

void free_sol( void *data, void *user_data )
{
    solution *sol = ( solution * )data;
    solution_free( sol );
    CC_IFFREE( sol, solution );
    sol = ( solution * )user_data;
}

void assign_iter( void *data, void *user_data )
{
    solution *sol = ( solution * )data;
    int *iter = ( int * )user_data;
    sol->iter = *iter;
}
void print_sol( void *data, void *user_data )
{
    solution *sol = ( solution * )data;
    int *i = ( int * )user_data;
    printf( "solution %d with totalweightcomptime %d \n", ( *i )++, sol->totalweightcomptime );
    solution_print( sol );
}

void print_totalweightcomptime( void *data, void *user_data )
{
    solution *sol  = ( solution * )data;
    int *i = ( int * )user_data;
    printf( "solution %d with totalweightcomptime %d, distance %d \n", ( *i )++, sol->totalweightcomptime,
            sol->dist );
}

void refset_dist( void *data, void *user_data )
{
    SSrefset_distance( ( SS * )user_data, ( solution * )data );
}

void for_each_comp_fitness( void *data, void *user_data )
{
    solution *sol = ( solution * ) data;
    int *lowerbound = ( ( int * )user_data );
    sol->fitness = ( float )( *lowerbound ) / ( sol->totalweightcomptime - *lowerbound );
}

/*
Compare functions
 */

int order_totalweightcomptime( const void *a, const void *b, void *data )
{
    (void) data;
    const int *aa = &( ( ( const solution * )a )->totalweightcomptime );
    const int *bb = &( ( ( const solution * )b )->totalweightcomptime );
    return *aa - *bb;
}
int order_distance( const void *a, const void *b, void *data )
{
    (void) data;
    const int *aa = &( ( ( const solution * )a )->dist );
    const int *bb = &( ( ( const solution * )b )->dist );
    return ( *aa - *bb );
}

int order_totalweightcomptime_list( const void *a, const void *b )
{
    const int *aa = &( ( ( const solution * )a )->totalweightcomptime );
    const int *bb = &( ( ( const solution * )b )->totalweightcomptime );
    return *aa - *bb;
}

int sort_vertices( const void *a, const void *b, void *data )
{
    (void) data;
    const int *aa = &( ( ( const Job * )a )->weight );
    const int *bb = &( ( ( const Job * )b )->weight );
    return -( *aa - *bb );
}

/*
Init and free functions
 */
void SS_init( SS *problem, int b1, int b2, double timelimit )
{
    problem->p         = ( P * )NULL;
    problem->rs        = ( REFSET * )NULL;
    problem->joblist = ( Job ** )NULL;
    problem->b1        = b1;
    problem->b2        = b2;
    problem->timelimit = timelimit;
    problem->nmachines = 0;
    problem->njobs     = 0;
    problem->status    = init;
    problem->random    = ( GRand * ) NULL;
    problem->iter = 0;
    problem->combine_method = 0;
}

void SS_free( SS *problem )
{
    REFSET_free( problem->rs );
    P_free( problem->p );
    CC_IFFREE( problem->rs, REFSET );
    CC_IFFREE( problem->p, P );
    problem->b2        = 0;
    problem->b1        = 0;
    problem->timelimit = .0;
    problem->njobs        = 0;

    if ( problem->random != NULL ) {
        g_rand_free( problem->random );
    }
}

void REFSET_init( REFSET *rs )
{
    if ( rs ) {
        rs->newsol        = 1;
        rs->list1         = g_queue_new();
        rs->list2         = g_queue_new();
    }
}

void REFSET_free( REFSET *rs )
{
    if ( rs ) {
        rs->newsol = 0;
        g_queue_foreach( rs->list1, free_sol, NULL );
        g_queue_foreach( rs->list2, free_sol, NULL );
        g_queue_free( rs->list1 );
        g_queue_free( rs->list2 );
        rs->newsol = 1;
        rs->list1 = ( GQueue * ) NULL;
        rs->list2 = ( GQueue * ) NULL;
    }
}

void P_init( P *p )
{
    if ( p ) {
        p->PSize = 0;
        p->list = ( GList * )NULL;
    }
}

void P_free( P *p )
{
    if ( p ) {
        for ( GList *it = p->list; it; it = it->next ) {
            solution_free( ( solution * )it->data );
            CC_IFFREE( it->data, solution );
        }

        p->PSize = 0;
        g_list_free( p->list );
        P_init( p );
    }
}

void free_list2( REFSET *rs )
{
    if ( rs ) {
        g_queue_foreach( rs->list2, free_sol, NULL );
        g_queue_free( rs->list2 );
        rs->list2 = ( GQueue * )NULL;
    }
}

/*
Scatter Search functions
 */

int SSproblem_definition(
    SS *problem,
    int b1,
    int b2,
    double timelimit,
    int combine_method,
    int njobs,
    int nmachines,
    Job *joblist,
    int lowerbound)
{
    int i, val         = 0;
    REFSET *temp_rs = NULL;
    P *temp_p       = NULL;

    /*initialize scatter search data structure */
    SS_init( problem, b1, b2, timelimit );
    problem->combine_method = combine_method;
    problem->nmachines = nmachines;
    problem->njobs = njobs;
    problem->lowerbound = lowerbound;
    problem->upperbound = INT_MAX;

    /* Initialize pool */
    problem->p = CC_SAFE_MALLOC(1, P);
    CCcheck_NULL_2( problem->p, "Failed to allocate memory to problem->p" );
    temp_p  = problem->p;
    P_init( temp_p );

    /* Initialize refset */
    problem->rs         = CC_SAFE_MALLOC( 1, REFSET );
    CCcheck_NULL_2( problem->rs, "Failed to allocate memory to problem->rs" );
    temp_rs             = problem->rs;
    REFSET_init( temp_rs );

    /* Initialize Joblist of scatter search data structure */
    problem->joblist = CC_SAFE_MALLOC(njobs, Job * );
    CCcheck_NULL_2( problem->joblist, "Failed to allocate memory" );
    for ( i = 0; i < njobs; i++ ) {
        problem->joblist[i] = joblist + i ;
    }

    problem->random = g_rand_new_with_seed( 48654642 );
    CCcheck_NULL_2( problem->random, "Failed in g_rand_new_with_seed" );

    CLEAN:
    if ( val ) {
        SS_free( problem );
    }
    return val;
}


void add_solution_pool( SS *scatter_search, solution *new_sol )
{
    P *p          = scatter_search->p;
    p->list = g_list_append( p->list, new_sol );
}

int add_solution_refset( SS *scatter_search)
{
    int val        = 0;
    REFSET *refset = scatter_search->rs;
    P *pool        = scatter_search->p;

    switch ( scatter_search->status ) {
        case init:
            pool->list = g_list_sort( pool->list, order_totalweightcomptime_list );
            for ( int i = 0; i < scatter_search->b2; ++i ) {
                void *data = pool->list->data;
                pool->list = g_list_remove( pool->list, data );
                g_queue_push_tail( refset->list1, data );
            }

            scatter_search->status = add;
            break;

        case add:
            for ( int i = 0; i < scatter_search->b1; ++i ) {
                void *data = maximum_distance( scatter_search );
                CCcheck_NULL_2( data, "Failed in maximum_distance" );
                pool->list = g_list_remove( pool->list, data );
                update_distance( scatter_search, ( solution * )data );
                g_queue_push_head( refset->list2, data );
            }

            scatter_search->status = update;
            break;

        case update:
            diversification_update( scatter_search );
            break;

        case opt:
            break;
    }

CLEAN:
    return val;
}

int solution_in_pool( SS *scatter_search, solution *new_sol )
{
    int val = 1;
    GList *it = scatter_search->p->list;
    int njobs = new_sol->njobs;

    if ( it == ( GList * ) NULL || g_list_length( it ) == 0 ) {
        val = 0;
        return val;
    }

    do {
        int i = 0;

        for ( i = 0; i < njobs; ++i ) {
            if ( ( ( solution * )it->data )->perm[i] != new_sol->perm[i] ) {
                break;
            }
        }

        if ( i == njobs ) {
            return val;
        }

        it = it->next;
    } while ( it != ( GList * ) NULL );

    val = 0;
    return val;
}

int solution_in_refset( SS *scatter_search, solution *new_sol )
{
    int val    = 1;
    GList *it  = scatter_search->rs->list1->head;
    int njobs = new_sol->njobs;
    int i = 0;

    if ( it == ( GList * ) NULL ) {
        val = 0;
        return val;
    }

    do {
        for ( i = 0; i < njobs; ++i ) {
            if ( ( ( solution * )it->data )->perm[i] != new_sol->perm[i] ) {
                break;
            }
        }

        if ( i == njobs ) {
            return val;
        }

        it = it->next;
    } while ( it != ( GList * ) NULL );

    it = scatter_search->rs->list2->head;

    if ( it == ( GList * ) NULL ) {
        val = 0;
        return val;
    }

    do {
        for ( i = 0; i < njobs; ++i ) {
            if ( ( ( solution * )it->data )->perm[i] != new_sol->perm[i] ) {
                break;
            }
        }

        if ( i == njobs ) {
            return val;
        }

        it = it->next;
    } while ( it != ( GList * ) NULL );

    val = 0;
    return val;
}

void print_pool( SS *scatter_search )
{
    P *pool = scatter_search->p;
    int i = 0;
    g_list_foreach( pool->list, print_sol, &i );
}

void print_pool_totalweightcomptime( SS *scatter_search )
{
    P *pool = scatter_search->p;
    int i = 0;
    g_list_foreach( pool->list, print_totalweightcomptime, &i );
}


void print_refset( SS *scatter_search )
{
    REFSET *rs = scatter_search->rs;
    int i = 0;
    g_queue_foreach( rs->list1, print_sol, &i );
    i = 0;
    g_queue_foreach( rs->list2, print_sol, &i );
}

void print_refset_totalweightcomptime( SS *scatter_search )
{
    REFSET *rs = scatter_search->rs;
    int i = 0;
    g_queue_foreach( rs->list1, print_totalweightcomptime, &i );
    g_queue_foreach( rs->list2, print_totalweightcomptime, &i );
}

void print_pool_n( SS *scatter_search, int n )
{
    P *pool = scatter_search->p;
    solution *sol = ( solution * )NULL;

    if ( ( guint )n > g_list_length( pool->list ) ) {
        printf( "n is too big\n" );
    } else {
        for ( int i = 0; i < n; ++i ) {
            sol = ( ( solution * )g_list_nth( pool->list, i )->data );
            printf( "solution %d with totalweightcomptime %d\n", i, sol->totalweightcomptime );
            solution_print( sol );
        }
    }
}

void print_list1( SS *scatter_search )
{
    REFSET *refset = scatter_search->rs;
    int i = 0;
    g_queue_foreach( refset->list1, print_sol, &i );
}

void print_distance( SS *scatter_search )
{
    P *pool = scatter_search->p;
    GList *it = ( GList * )NULL;
    int i = 0;

    for ( it = pool->list; it; it = it->next ) {
        printf( "Distance to refset = %d and solution %d\n",
                ( ( solution * )it->data )->dist, i++ );
    }
}


int SSrefset_distance( SS *scatter_search, solution *new_sol )
{
    int val = 0;
    REFSET *refset = scatter_search->rs;
    min_max temp = {scatter_search->njobs, INT_MAX, new_sol};

    if ( refset->list1->head == ( GList * )NULL ) {
        printf( "We can't compute the distance between pool and refset\n" );
        val = 1;
        return val;
    }

    g_queue_foreach( refset->list1, distance_min_max, &temp );
    g_queue_foreach( refset->list2, distance_min_max, &temp );
    new_sol->dist = temp.min;
    
    return val;
}



void *maximum_distance( SS *scatter_search )
{
    solution *val = ( solution * )NULL;
    P *pool = scatter_search->p;
    val = ( solution * )pool->list->data;
    g_list_foreach( pool->list, max_dist, val );
    return val;
}

int update_distance( SS *scatter_search, solution *sol )
{
    int l, k, val = 0;
    GList *it = ( GList * )NULL;
    P *pool = scatter_search->p;

    if ( pool->list == ( GList * )NULL ) {
        printf( "We can't update the distances. The pool is empty!!!!\n" );
        val = 1;
        goto CLEAN;
    }

    for ( it = pool->list; it; it = g_list_next( it ) ) {
        int curval = 0;

        for ( l = 0; l < scatter_search->njobs - 1; l++ ) {
            for ( k = l + 1; k < scatter_search->njobs; k++ ) {
                if ( ( ( ( solution * )it->data )->vlist[l].part == ( ( solution * )
                        it->data )->vlist[k].part && sol->vlist[l].part != sol->vlist[k].part ) ||
                        ( ( ( solution * )it->data )->vlist[l].part != ( ( solution * )
                                it->data )->vlist[k].part && sol->vlist[l].part == sol->vlist[k].part ) ) {
                    curval++;
                }

                if ( curval > ( ( solution * )it->data )->dist ) {
                    break;
                }
            }

            if ( curval > ( ( solution * )it->data )->dist ) {
                break;
            }
        }

        if ( curval <= ( ( solution * )it->data )->dist ) {
            ( ( solution * )it->data )->dist = curval;
        }
    }

CLEAN:
    return val;
}



int SSCreate_refset( SS *scatter_search )
{
    int val   = 0;
    P *pool   = scatter_search->p;
    add_solution_refset( scatter_search );
    g_list_foreach( pool->list, refset_dist, scatter_search );
    add_solution_refset( scatter_search );

    return val;
}

int compute_fitness( SS *scatter_search )
{
    int val = 0;
    REFSET *refset = scatter_search->rs;
    int lowerbound = scatter_search->lowerbound;
    g_queue_foreach( refset->list1, for_each_comp_fitness, &lowerbound );
    g_queue_foreach( refset->list2, for_each_comp_fitness, &lowerbound );
    return val;
}

int SSrun_scatter_search( SS *scatter_search)
{
    int flag, val         = 0;
    REFSET *refset    = scatter_search->rs;
    int nbnew_sol = 0;
    int nb_noimprovements = 0;
    int *totsubset = ( int * ) NULL;
    totsubset = CC_SAFE_MALLOC( scatter_search->b2 + 1, int );
    CCcheck_NULL_2( totsubset, "Failed tot allocate memory to tot subset" );

    if ( refset->list1->head == ( GList * )NULL ) {
        printf( "We can't run scatter search, refset is empty\n" );
        val = 1;
        goto CLEAN;
    }

    while ( refset->newsol && scatter_search->status != opt) {
        GQueue *list = g_queue_copy( refset->list1 );

        for ( GList *it = refset->list2->head; it; it = it->next ) {
            g_queue_push_tail( list, it->data );
        }

        compute_fitness( scatter_search );
        refset->newsol = 0;

        if ( scatter_search->iter >= 0 ) {
            int best = ( ( solution * )refset->list1->head->data )->totalweightcomptime;
            printf( "iteration %d with best totalweightcomptime %d, number of new solutions %d\n",
                    scatter_search->iter, best, nbnew_sol );

            if ( nbnew_sol == 0 ) {
                nb_noimprovements++;
            } else {
                nb_noimprovements = 0;
            }

            nbnew_sol = 0;
        }

        k_subset_init( g_queue_get_length( list ), 2, totsubset, &flag );

        while ( flag && scatter_search->status != opt ) {
            solution new_sol;
            solution_init( &new_sol );
            int rval = 1;

            if ( scatter_search->combine_method == 0 ) {
                rval = combine_GPX( scatter_search, list, totsubset, 2, &new_sol );
            } else {
                rval = combine_PM( scatter_search, list, totsubset, 2, &new_sol );
            }

            if ( !rval ) {
                /* check feasibility schedule */
                localsearch_random_k( &new_sol, scatter_search->lowerbound,2 );
                solution_unique( &new_sol );

                if ( !dynamic_update( scatter_search, list, &new_sol ) ) {
                    if ( new_sol.totalweightcomptime < scatter_search->upperbound ) {
                        scatter_search->upperbound = new_sol.totalweightcomptime;
                    }

                    if ( new_sol.totalweightcomptime == scatter_search->lowerbound ) {
                        scatter_search->status = opt;
                        printf( "Found optimal with SS\n" );
                    }
                    nbnew_sol++;
                }
            }

            solution_free( &new_sol );
            k_subset_lex_successor( g_queue_get_length( list ), 2, totsubset, &flag );
        }

        k_subset_init( g_queue_get_length( list ), 3, totsubset, &flag );
        solution *sol = ( solution * )NULL;
        sol = g_queue_peek_nth( list, totsubset[1] - 1 );

        while ( flag && scatter_search->status != opt
                && sol->totalweightcomptime <= scatter_search->upperbound ) {
            solution new_sol;
            solution_init( &new_sol );
            int rval = 1;

            if ( scatter_search->combine_method == 0 ) {
                rval = combine_GPX( scatter_search, list, totsubset, 3, &new_sol ); //GPX
            } else {
                rval = combine_PM( scatter_search, list, totsubset, 3,
                                   &new_sol ); //Dell'Amico et al.
            }

            if ( !rval ) {
                /* Check solution */
                CCcheck_val( val, "Failed in solution_coloring" );
                localsearch_random_k( &new_sol, scatter_search->lowerbound,3 );
                solution_unique( &new_sol );

                if ( !dynamic_update( scatter_search, list, &new_sol ) ) {
                    if ( new_sol.totalweightcomptime < scatter_search->upperbound ) {
                        scatter_search->upperbound = new_sol.totalweightcomptime;
                    }

                    if ( new_sol.totalweightcomptime == scatter_search->lowerbound ) {
                        scatter_search->status = opt;
                        printf( "Found optimal with SS\n" );
                    }

                    nbnew_sol++;
                }
            }

            solution_free( &new_sol );
            k_subset_lex_successor( g_queue_get_length( list ), 3, totsubset, &flag );
            sol = g_queue_peek_nth( list, totsubset[1] - 1 );
        }

        k_subset_init( g_queue_get_length( list ), 4, totsubset, &flag );
        sol = g_queue_peek_nth( list, totsubset[1] - 1 );
        solution *temp_sol = g_queue_peek_nth( list, totsubset[2] - 1 );

        while ( flag && scatter_search->status != opt
                && sol->totalweightcomptime <= scatter_search->upperbound
                && temp_sol->totalweightcomptime <= scatter_search->upperbound ) {
            solution new_sol;
            solution_init( &new_sol );
            int rval = 1;

            if ( scatter_search->combine_method == 0 ) {
                rval = combine_GPX( scatter_search, list, totsubset, 4, &new_sol );
            } else {
                rval = combine_PM( scatter_search, list, totsubset, 4, &new_sol );
            }

            if ( !rval ) {
                /* Check solution */
                CCcheck_val( val, "Failed in solution_coloring" );
                localsearch_random_k(&new_sol, scatter_search->lowerbound,3 );
                solution_unique( &new_sol );

                if ( !dynamic_update( scatter_search, list, &new_sol ) ) {
                    if ( new_sol.totalweightcomptime < scatter_search->upperbound ) {
                        scatter_search->upperbound = new_sol.totalweightcomptime;
                    }

                    if ( new_sol.totalweightcomptime == scatter_search->lowerbound ) {
                        scatter_search->status = opt;
                        printf( "Found optimal with SS\n" );
                    }

                    nbnew_sol++;
                }
            }

            solution_free( &new_sol );
            k_subset_lex_successor( g_queue_get_length( list ), 4, totsubset, &flag );
            sol = g_queue_peek_nth( list, totsubset[1] - 1 );
            temp_sol = g_queue_peek_nth( list, totsubset[2] - 1 );
        }

        for ( int i = 0; i < scatter_search->b2 + 1; ++i ) {
            totsubset[i] = i;
        }

        for ( int i = 5; i < scatter_search->b2 + 1
                && scatter_search->status != opt; i++ ) {
            solution new_sol;
            solution_init( &new_sol );
            int rval = 1;

            if ( scatter_search->combine_method == 0 ) {
                rval = combine_GPX( scatter_search, list, totsubset, i, &new_sol );
            } else {
                rval = combine_PM( scatter_search, list, totsubset, i, &new_sol );
            }

            if ( !rval ) {
                localsearch_random_k(&new_sol, scatter_search->lowerbound, 3 );
                solution_unique( &new_sol );

                if ( !dynamic_update( scatter_search, list, &new_sol ) ) {
                    if ( new_sol.totalweightcomptime < scatter_search->upperbound ) {
                        scatter_search->upperbound = new_sol.totalweightcomptime;
                    }

                    if ( new_sol.totalweightcomptime == scatter_search->lowerbound ) {
                        scatter_search->status = opt;
                        printf( "Found optimal with SS\n" );
                    }

                    nbnew_sol++;
                }

                solution_free( &new_sol );
            }
        }

        g_queue_free( list );
        scatter_search->iter++;

        if ( refset->newsol == 0 && scatter_search->iter < 2 * scatter_search->njobs
                 && nb_noimprovements < 10 ) {
            add_solution_refset( scatter_search );
            refset->newsol = 1;
        }

    }

CLEAN:
    CC_IFFREE( totsubset, int );
    return val;
}

/* functions for constructing new solutions*/

static int select_parrent( GRand *Rand, int nbelements, int *sum, int totsum )
{
    int val = -1;
    int counter = 0;
    int *parents = ( int * )NULL;

    if ( totsum <= 0 ) {
        goto CLEAN;
    }

    for ( int i = 0; i < nbelements; ++i ) {
        if ( sum[i] > 0 ) {
            counter++;
        }
    }

    if ( counter == 0 ) {
        goto CLEAN;
    }

    int temp = counter - 1;
    parents = CC_SAFE_MALLOC( counter, int );

    for ( int i = 0; i < nbelements; ++i ) {
        if ( sum[i] > 0 ) {
            parents[temp--] = i;
        }
    }

    val = parents[g_rand_int_range( Rand, 0, counter )];
CLEAN:
    CC_IFFREE( parents, int );
    return val;
}

int combine_GPX( SS *scatter_search, GQueue *queue, int *subsetsol,
                 int nbelements, solution *new_sol )
{
    int k, val    = 1;
    int nmachines = ( ( solution * )scatter_search->rs->list1->head->data )->nmachines;
    int njobs = ( ( solution * )scatter_search->rs->list1->head->data )->njobs;
    int part, i, j;
    int totsum = 0;
    int *sum = ( int * ) NULL;
    solution *sol = ( solution * )NULL;
    solution *temp_sol = ( solution * ) NULL;
    sum = CC_SAFE_MALLOC( nbelements, int );
    CCcheck_NULL_2( sum, "Failed to allocate memory" );

    for ( i = 0; i < nbelements; i++ ) {
        sum[i] = njobs;
        totsum += njobs;
    }

    for ( k = 1;  k <= nbelements && val; k++ ) {
        temp_sol  = g_queue_peek_nth( queue, subsetsol[k] - 1 );

        if ( temp_sol->iter >= scatter_search->iter ) {
            val = 0;
        }
    }

    if ( val ) {
        goto CLEAN;
    }

    sol = CC_SAFE_MALLOC( nbelements, solution );
    CCcheck_NULL_2( sol, "Failed to allocate memory" );

    for ( i = 0; i  < nbelements; i++ ) {
        solution_copy( sol + i, *( ( solution * )g_queue_peek_nth( queue,
                                   subsetsol[i + 1] - 1 ) ) );
    }

    val = solution_alloc( new_sol, nmachines,  njobs );
    CCcheck_val_2( val, "Failed in solution_alloc" );
    new_sol->nmachines = 0;
    new_sol->njobs = 0;

    while ( new_sol->nmachines < nmachines && totsum > 0 ) {
        i = select_parrent( scatter_search->random, nbelements, sum, totsum );

        if ( i != -1 ) {
            part = minimum_partition( sol + i, scatter_search->lowerbound );

            for ( GList *it = ( sol + i )->part[part].list->head; it; it = it->next ) {
                partlist_insert( &new_sol->part[new_sol->nmachines], new_sol->vlist,
                                 ( Job * )it->data );
            }

            for ( GList *it = new_sol->part[new_sol->nmachines].list->head; it ;
                    it = it->next ) {
                for ( j = 0; j < nbelements; j++ ) {
                    partlist_delete_custom( ( sol + j )->vlist, ( ( Job * )it->data ), njobs );
                }
            }

            for ( j = 0; j < nbelements; j++ ) {
                sum[j] -= g_queue_get_length( new_sol->part[new_sol->nmachines].list );
            }

            totsum -= nbelements * g_queue_get_length(
                          new_sol->part[new_sol->nmachines].list );
            new_sol->njobs += g_queue_get_length( new_sol->part[new_sol->nmachines].list );
            new_sol->nmachines++;
        } else {
            solution_print( new_sol );

            for ( j = 0; j < nbelements; j++ ) {
                printf( "%d ", sum[j] );
            }

            printf( "\n" );
            printf( "totsum = %d\n", totsum );
        }
    }

    if ( new_sol->njobs != njobs ) {
        pmcheap *heap = ( pmcheap * )NULL;
        partlist *temp = (partlist*) NULL;
        GQueue *list = g_queue_new();

        for ( i = 0; i < njobs; ++i ) {
            if ( sol->part[i].list != NULL ) {
                for ( GList *it = sol->part[i].list->head; it; it = it->next ) {
                    g_queue_push_head( list, it->data );
                }
            }
        }

        g_queue_sort( list, sort_vertices, NULL );
        pmcheap_init( &heap, new_sol->nmachines );

        for ( i = 0; i < new_sol->nmachines; ++i ) {
            val = pmcheap_insert( heap, new_sol->part[i].completiontime, new_sol->part + i );

            if ( val ) {
                printf( "Failed at pmcheap_insert in %s at %d\n", __FILE__, __LINE__ );
                g_queue_free( list );
                pmcheap_free( heap );
                goto CLEAN;
            }
        }

        GList *it = list->head;

        while ( it ) {
            temp = pmcheap_min( heap );
            Job *job = ( Job * )it->data;


            partlist_insert( temp, new_sol->vlist, job);
            new_sol->njobs++;
            pmcheap_insert( heap, temp->completiontime, temp );

            it = it->next;
        }

        pmcheap_free( heap );
        g_queue_free( list );
    }

CLEAN:

    if ( val ) {
        solution_free( new_sol );
    }

    CC_IFFREE( sum, int );

    if ( sol ) {
        for ( i = 0; i < nbelements; i++ ) {
            solution_free( sol + i );
        }

        CC_IFFREE( sol, solution );
    }

    return val;
}

int combine_PM( SS *scatter_search, GQueue *queue, int *subsetsol,
                int nbelements, solution *new_sol )
{
    int i, j, k, winner, val = 1;
    int nbsubsets = ( ( scatter_search->njobs  + 1 ) * scatter_search->njobs ) / 2;
    int njobs, nmachines, it, count, step, first;
    pmcheap *heap = ( pmcheap * ) NULL;
    pmcheap *nodeheap = ( pmcheap * )NULL;
    partlist **temp = ( partlist ** )NULL;
    solution *sol = ( solution * )NULL;
    float *fitness = ( float * )NULL;
    float *cumulfitness = ( float * ) NULL;

    if ( g_queue_get_length( queue ) >= ( guint )nbelements ) {
        njobs = ( ( solution * )queue->head->data )->njobs;
        nmachines = ( ( solution * )queue->head->data )->nmachines;
    } else {
        printf( "Error in combine solution.\n" );
        val = 1;
        goto CLEAN;
    }

    for ( k = 1;  k <= nbelements && val; k++ ) {
        sol  = g_queue_peek_nth( queue, subsetsol[k] - 1 );

        if ( sol->iter >= scatter_search->iter ) {
            val = 0;
        }
    }

    if ( val ) {
        goto CLEAN;
    }

    temp = CC_SAFE_MALLOC( nmachines, partlist * );
    CCcheck_NULL_2( temp, "Failed to allocate memory" );
    fitness = CC_SAFE_MALLOC( nbsubsets, float );
    CCcheck_NULL_2( fitness, "Failed to allocate memory" );
    fill_float( fitness, nbsubsets, .0 );
    cumulfitness = CC_SAFE_MALLOC( nbsubsets, float );
    CCcheck_NULL_2( cumulfitness, "Failed to allocate memory" );
    fill_float( cumulfitness, nbsubsets, .0 );
    val = solution_alloc( new_sol, nmachines, njobs );
    CCcheck_val_2( val, "Failed in solution_alloc" );
    new_sol->njobs = 0;
    val = pmcheap_init( &heap, nmachines );
    CCcheck_val_2( val, "Failed in pmcheap_init" );

    for ( i = 0; i < nmachines; i++ ) {
        val = pmcheap_insert( heap, new_sol->part[i].completiontime, new_sol->part + i );
        CCcheck_val_2( val, "Failed in pmcheap_insert" );
    }

    val = pmcheap_init( &nodeheap, njobs );
    CCcheck_NULL_2( nodeheap, "Failed in pmcheap_init" );

    for ( i = 0; i < scatter_search->njobs; i++ ) {
        for ( j = i ; j < scatter_search->njobs; j++ ) {
            int key = nodepair_ref_key( i, j );
            fitness[key] = 0;

            if ( i != j ) {
                for ( k = 1; k <= nbelements; k++ ) {
                    sol  = g_queue_peek_nth( queue, subsetsol[k] - 1 );

                    if ( sol->vlist[i].part == sol->vlist[j].part ) {
                        fitness[key] += sol->fitness;
                    }
                }
            }
        }
    }

    for ( i = 1; i < nbsubsets; i++ ) {
        cumulfitness[i] = fitness[i] + cumulfitness[i - 1];
    }

    while ( cumulfitness[nbsubsets - 1] > 0 ) {
        const float f = g_rand_double_range( scatter_search->random, .0,
                                             cumulfitness[nbsubsets - 1] );
        count = nbsubsets;
        first = 0;

        while ( count > 0 ) {
            it = first;
            step = count / 2;
            it += step;

            if ( cumulfitness[it] < f ) {
                first = ++it;
                count -= step + 1;
            } else {
                count = step;
            }
        }

        winner = first;
        inodepair_ref_key( &i, &j, winner );
        int counter = 0;
        temp[counter] = pmcheap_min( heap );
        Job *node1 = *( scatter_search->joblist + i );
        Job *node2 = *( scatter_search->joblist + j );


    
        val = partlist_insert( temp[counter], new_sol->vlist, node1 );
        val = partlist_insert( temp[counter], new_sol->vlist, node2 );
        CCcheck_val_2( val, "Failed partlist_insert" );
        CCcheck_val_2( val, "Failed partlist_insert" );
        new_sol->njobs += 2;


        if ( counter >= nmachines ) {
            counter--;
        }

        for ( i = 0; i <= counter; i++ ) {
            pmcheap_insert( heap, temp[i]->completiontime, temp[i] );
        }

        for ( i = 0; i < scatter_search->njobs; i++ ) {
            if ( i <= node1->job ) {
                fitness[nodepair_ref_key( i, node1->job )] = .0;
            } else {
                fitness[nodepair_ref_key( node1->job, i )] = .0;
            }

            if ( i <= node2->job ) {
                fitness[nodepair_ref_key( i, node2->job )] = .0;
            } else {
                fitness[nodepair_ref_key( node2->job, i )] = .0;
            }
        }

        for ( i = 1; i < nbsubsets; i++ ) {
            cumulfitness[i] = fitness[i] + cumulfitness[i - 1];
        }
    }

    if ( new_sol->njobs != njobs ) {
        for ( i = 0; i < njobs; i++ ) {
            if ( new_sol->vlist[scatter_search->joblist[i]->job].part == NULL ) {
                pmcheap_insert( nodeheap, -( scatter_search->joblist[i]->processingime ),
                                scatter_search->joblist[i] );
            }
        }

        pmcheap_reset( heap );

        for ( i = 0; i < nmachines; i++ ) {
            pmcheap_insert( heap, new_sol->part[i].completiontime, new_sol->part + i );
        }

        while ( nodeheap->end > 0 ) {
            Job *node = pmcheap_min( nodeheap );
            CCcheck_NULL_2( node, "Failed pmcheap_min" );
            int counter = 0;
            temp[counter] = pmcheap_min( heap );


            partlist_insert( temp[counter], new_sol->vlist, node );
            new_sol->njobs++;

            if ( counter >= nmachines ) {
                counter--;
            }

            for ( i = 0; i <= counter; i++ ) {
                pmcheap_insert( heap, temp[i]->completiontime, temp[i] );
            }
        }
    }

CLEAN:

    if ( val ) {
        solution_free( new_sol );
    }

    pmcheap_free( nodeheap );
    pmcheap_free( heap );
    CC_IFFREE( temp, partlist * );
    CC_IFFREE( fitness, float );
    CC_IFFREE( cumulfitness, float );
    return val;
}



int static_update( SS *scatter_search )
{
    int val = 0;
    GList *it = ( GList * )NULL;
    REFSET *refset = scatter_search->rs;
    P *pool = scatter_search->p;

    for ( it = refset->list1->head; it; it = it->next ) {
        ( ( solution * )it->data )->iter = 1;
    }

    for ( it = refset->list2->head; it; it = it->next ) {
        ( ( solution * )it->data )->iter = 1;
    }

    for ( it = pool->list; it; it = it->next ) {
        SSrefset_distance( scatter_search, it->data );
    }

    refset->newsol = 0;
    it = pool->list;

    while ( it ) {
        solution *sol = ( solution * )it->data;
        solution *last = g_queue_peek_tail( refset->list1 );
        solution *last1 = g_queue_peek_tail( refset->list2 );
        int not_in_refset = !solution_in_refset( scatter_search, sol );

        if ( sol->totalweightcomptime < last->totalweightcomptime && not_in_refset ) {
            refset->newsol = 1;
            void *data = g_queue_pop_tail( refset->list1 );
            solution_free( data );
            CC_IFFREE( data, solution );
            it = pool->list = g_list_remove( pool->list, sol );

            if ( g_list_length( pool->list ) != 0 ) {
                val = update_distance( scatter_search, sol );
                CCcheck_val_2( val, "Failed in update_distance" );
                g_queue_insert_sorted( refset->list1, sol, order_totalweightcomptime, NULL );
            }
        } else
            if ( sol->dist > last1->dist ) {
                refset->newsol = 1;
                solution *data = g_queue_pop_tail( refset->list2 );
                solution_free( data );
                CC_IFFREE( data, solution );
                it = pool->list = g_list_remove( pool->list, sol );

                if ( g_list_length( pool->list ) != 0 ) {
                    val = update_distance( scatter_search, sol );
                    CCcheck_val_2( val, "Failed in update_distance" );
                    g_queue_insert_sorted( refset->list2, sol, order_distance, NULL );
                }
            } else {
                it = it->next;
            }
    }

    printf( "update\n" );
CLEAN:
    return val;
}

int dynamic_update( SS *scatter_search, GQueue *list, solution *new_sol )
{
    int val = 1;
    REFSET *refset = scatter_search->rs;
    solution *last1 = g_queue_peek_tail( refset->list1 );
    solution *last2 = g_queue_peek_head( refset->list2 );
    int not_in_refset = !solution_in_refset( scatter_search, new_sol );

    if ( new_sol->totalweightcomptime <  last1->totalweightcomptime && not_in_refset ) {
        refset->newsol = 1;
        solution *data = ( solution * )g_queue_pop_tail( refset->list1 );
        solution *data1 = ( solution * )g_queue_find( list, data )->data;
        int n = g_queue_index( list, data1 );
        g_queue_pop_nth( list, n );
        solution_update( data, *( new_sol ) );
        data->iter = scatter_search->iter + 1;
        g_queue_insert_sorted( refset->list1, data, order_totalweightcomptime, NULL );
        g_queue_insert_sorted( list, data, order_totalweightcomptime, NULL );
        val = 0;
        return val;
    } else {
        SSrefset_distance( scatter_search, new_sol );

        if ( new_sol->dist > last2->dist ) {
            refset->newsol = 1;
            solution *data = g_queue_pop_head( refset->list2 );
            solution *data1 = ( solution * )g_queue_find( list, data )->data;
            int n = g_queue_index( list, data1 );
            g_queue_pop_nth( list, n );
            solution_update( data, *( new_sol ) );
            data->iter = scatter_search->iter + 1;
            g_queue_insert_sorted( refset->list2, data, order_distance, NULL );
            g_queue_insert_sorted( list, data, order_distance, NULL );
            val = 0;
            return val;
        }
    }

    return val;
}

int diversification_update( SS *scatter_search )
{
    int val = 0;
    REFSET *refset = scatter_search->rs;
    free_list2( refset );
    refset->list2 = g_queue_new();

    while ( ( int )g_queue_get_length( refset->list2 ) < scatter_search->b1
            && scatter_search->status != opt  ) {

        if ( g_rand_int_range( scatter_search->random, 0, 2 ) == 0 ) {
            if ( g_rand_int_range( scatter_search->random, 0, 2 ) == 0 ) {
                    /* Construct feasible schedules */
            } else {
               /* Construct feasible schedules */
            }

            if ( val /* Check if new solution */) {
                solution *new_sol = ( solution * )NULL;
                new_sol = CC_SAFE_MALLOC( 1, solution );
                CCcheck_NULL( new_sol, "Failed to allocate memory to new_sol" );
                solution_init( new_sol );
                solution_alloc( new_sol, scatter_search->nmachines, scatter_search->njobs );
                localsearch_SS( new_sol, scatter_search->lowerbound );
                solution_unique( new_sol );
                /* Check feasibility */
                CCcheck_val( val, "Failed in solution_coloring" );

                if ( !solution_in_refset( scatter_search, new_sol ) ) {
                    if ( new_sol->totalweightcomptime < scatter_search->upperbound ) {
                        scatter_search->upperbound = new_sol->totalweightcomptime;
                    }

                    if ( scatter_search->upperbound == scatter_search->lowerbound ) {
                        scatter_search->status = opt;
                    }

                    if ( new_sol->totalweightcomptime < ( ( solution * )g_queue_peek_tail(
                                                   refset->list1 ) )->totalweightcomptime ) {
                        void *data = g_queue_pop_tail( refset->list1 );
                        g_queue_insert_sorted( refset->list1, new_sol, order_totalweightcomptime, NULL );
                        SSrefset_distance( scatter_search, data );
                        g_queue_push_head( refset->list2, data );
                    } else {
                        SSrefset_distance( scatter_search, new_sol );
                        g_queue_push_head( refset->list2, new_sol );
                    }
                } else {
                    solution_free( new_sol );
                    CC_IFFREE( new_sol, solution );
                }
            }
        }
    }

    if ( g_queue_get_length( refset->list2 ) == 0 ) {
        printf( "NO new solutions\n" );
        val = 1;
        goto CLEAN;
    }

    g_queue_sort( refset->list2, order_distance, NULL );
    g_queue_foreach( refset->list1, assign_iter, &scatter_search->iter );
    g_queue_foreach( refset->list2, assign_iter, &scatter_search->iter );
CLEAN:
    return val;
}
