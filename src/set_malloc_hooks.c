/* Test program to measure the heap memory available for user allocation. */
/* Prototypes for __malloc_hook, __free_hook */
#include <malloc.h>
     
/* Prototypes for our hooks.  */
static void my_init_hook (void);
static void *my_malloc_hook (size_t, const void *);
static void my_free_hook (void*, const void *);
     
static void *(*old_malloc_hook)(size_t, const void *);
static void (*old_free_hook)(void*,const void*);

/* Override initializing hook from the C library. */
void (*__malloc_initialize_hook) (void) = my_init_hook;

#define MAX_MALLOC_BLOCKS 100000
const int debugging=0;
const int all_output=0;

typedef struct malloc_block_s {
  void *pointer;
  long int size;
} malloc_block_t;
/* blocks array contains the allocated blocks, index 0 is
   the address, index 1 the size. */
static malloc_block_t blocks[MAX_MALLOC_BLOCKS];
static long int last_block_index;
static long int actual_allocated;
static long int max_allocated,max_number_of_blocks;
static int just_count;

static void my_init_hook (void) {
  int i;
  char *value;

  /*printf("Initializing the malloc hooks.\n");*/
  old_malloc_hook = __malloc_hook;
  old_free_hook = __free_hook;
  /* read in environment variable MALLOC_HOOK_COUNT
     if set, only count the maximum number of simultanously
     set malloc blocks. */
  value = getenv("MALLOC_HOOK_COUNT");
  if (value==NULL) just_count=0;
  else {
    if (strtol(value,NULL,10)!=0) {
      just_count = 1;
    } else {
      just_count = 0;
    }
  }
  actual_allocated = 0L;
  last_block_index=-1;
  max_number_of_blocks=0;
  if (!just_count) {
    for (i=0;i<MAX_MALLOC_BLOCKS;i++) {
      blocks[i].pointer=NULL;
      blocks[i].size = 0L;
    }
  }
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
}
     
static void *my_malloc_hook (size_t size, const void *caller) {
  void *result;
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
  /* Call recursively */
  result = malloc (size);
  /* Save underlying hooks */
  old_malloc_hook = __malloc_hook;
  old_free_hook = __free_hook;
  /* printf might call malloc, so protect it too. */
  if (result!=NULL) {
    /* malloc was successful */
    actual_allocated += size;
    last_block_index++;
    if (last_block_index>max_number_of_blocks) 
      max_number_of_blocks=last_block_index;

    if (!just_count) {
      if (last_block_index<MAX_MALLOC_BLOCKS) {
	blocks[last_block_index].pointer=result;
	blocks[last_block_index].size=size;
      } else {
	printf("\n----------------\n last_block_index =%ld exceeds MAX_MALLOC_BLOCKS\n",last_block_index);
      }
      if (actual_allocated>max_allocated) {
	max_allocated=actual_allocated;
	/*printf("max_allocated = %lu\n",max_allocated);*/
      }
    }
  }
  /*printf("actual_allocate = %lu\n",actual_allocated);*/
  /*printf ("malloc (%u) returns %p\n", (unsigned int) size, result);*/
  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
  return result;
}

static void
my_free_hook (void *ptr, const void *caller)
{
  int i, act_index;
  long int size;
  /* Restore all old hooks */
  __malloc_hook = old_malloc_hook;
  __free_hook = old_free_hook;
  /* find the right entry in blocks, start with last_block_index
     and search backwards for the right ptr */
  if (ptr!=NULL) {
    if (!just_count) {
      size = -1;
      for (i=last_block_index;i>=0;i--) {
	if (blocks[i].pointer==ptr) {
	  size = blocks[i].size;
	  act_index = i;
	}
      }
    } else {
      /* set size to some value != -1 to call the free
	 function in the next block */
      size = 0;
    }
    if (size!=-1) {
      /* Call recursively */
      free(ptr);
      /* Save underlying hooks */
      old_malloc_hook = __malloc_hook;
      old_free_hook = __free_hook;
      
      if (!just_count) {
	actual_allocated -= size;
	/* copy all later blocks to overwrite the 
	   actual freed one. */
	for (i=act_index+1;i<=last_block_index;i++) {
	  blocks[i-1]=blocks[i];
	}
      }
      last_block_index--;
      /*printf("After freeing: actual_allocated = %lu\n",actual_allocated);*/
    } else {
      /* printf might call free, so protect it too. */
      if (all_output) {
	printf("Error: A block with address %p has never been allocated.\n",ptr);
	printf("actual_allocated = %lu, last_block_index = %ld\n",
	       actual_allocated,last_block_index);
	if (last_block_index>=0) {
	  for (i=0;i<=last_block_index;i++) {
	    if (blocks[i].pointer==NULL) break;
	    printf("(%p ==> %lu) ",blocks[i].pointer,blocks[i].size);
	  }
	  printf("\n");
	} else {
	  printf("last_block_index < 0, output blocks:\n");
	  for (i=0;i<MAX_MALLOC_BLOCKS;i++) {
	    if (blocks[i].pointer==NULL) break;
	    printf("(%p ==> %lu) ",blocks[i].pointer,blocks[i].size);
	  }
	  printf("\n");
	}
	printf("Max allocated is %lu.\n",max_allocated);
      }
    
      /*exit(1);*/
    }
  } else {
    if (debugging) printf("free has been called with ptr=NULL! Ignoring!\n");
  }

  /*printf ("freed pointer %p\n", ptr);*/
  /* Restore our own hooks */
  __malloc_hook = my_malloc_hook;
  __free_hook = my_free_hook;
}

void get_malloc_stat(long int *max_all,int *unfreed_blocks, long int *still_all) {
  *max_all=max_allocated;
  *unfreed_blocks = last_block_index+1;
  *still_all = actual_allocated;
}

void show_actual_blocks() {
  int i;
  for (i=0;i<=last_block_index;i++) {
    printf("%d. block, pointer = %p, size = %lu\n",i,blocks[i].pointer,blocks[i].size);
  }
}

long int get_allocated_memory(void) {
  return actual_allocated;
}

long int get_max_number_of_blocks(void) {
  return max_number_of_blocks;
}

