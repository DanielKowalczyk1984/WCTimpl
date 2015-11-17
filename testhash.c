//ex-ghashtable-1.c
#include <glib.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv) {
 GHashTable* hash = g_hash_table_new(g_direct_hash, g_direct_equal);
 g_hash_table_insert(hash, GINT_TO_POINTER(1), GINT_TO_POINTER(1));
 g_hash_table_insert(hash, GINT_TO_POINTER(2), GINT_TO_POINTER(1));
 g_hash_table_insert(hash, GINT_TO_POINTER(3), GINT_TO_POINTER(1));
 printf("There are %d keys in the hash\n", g_hash_table_size(hash));
 printf("The capital of Texas is %d\n", GPOINTER_TO_INT(g_hash_table_lookup(hash,GINT_TO_POINTER(2))));
 gboolean found = g_hash_table_remove(hash, "Virginia");
 printf("The value 'Virginia' was %sfound and removed\n", found ? "" : "not ");
 g_hash_table_destroy(hash);
 GList *list = (GList*) NULL;
 int i;
 for( i = 0; i < 10; ++i) {
     list = g_list_append(list, GINT_TO_POINTER(i));
 }
 GList *it;
 for(it = list; it; it = it->next) {
    int C = GPOINTER_TO_INT(it->data);
    int D = GPOINTER_TO_INT(it->next->data);
    if(C < D) {
        printf("Daniel\n");
    }
     printf("%d\n", C + D);
 }
 return 0;
}
