/**
 * @file ParallelHashMap.h
 * @brief A thread-safe hash map supporting insertion and lookup operations.
 * @details The parallel hash map is built on top of a fixed-sized hash map
 *    object and features OpenMP concurrency structures. The underlying
 *    fixed-sized hash map handles collisions with chaining.
 * @date June 6, 2015
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef __PARALLEL_HASH_MAP__
#define __PARALLEL_HASH_MAP__
#include<iostream>
#include<stdexcept>
#include<functional>
#include<omp.h>

#include "log.h"


/**
 * @class FixedHashMap ParallelHashMap.h "src/ParallelHashMap.h"
 * @brief A fixed-size hash map supporting insertion and lookup operations.
 * @details The FixedHashMap class supports insertion and lookup operations
 *    but not deletion as deletion is not needed in the OpenMOC application.
 *    This hash table uses chaining for collisions and does not incorporate
 *    concurrency objects except for tracking the number of entries in the
 *    table for which an atomic increment is used. This hash table is not
 *    thread safe but is used as a building block for the ParallelHashMap
 *    class. This table guarantees O(1) insertions and lookups on average.
 */
template <class K, class V>
class FixedHashMap {
  struct node {
    node(K k_in, V v_in) : next(NULL), key(k_in), value(v_in) {}
    K key;
    V value;
    node *next;
  };

  private:
    size_t _M;      /* table size */
    size_t _N;      /* number of elements present in table */
    node ** _buckets;   /* buckets of values stored in nodes */

  public:

    FixedHashMap(size_t M = 64);
    virtual ~FixedHashMap();
    bool contains(K& key);
    V& at(K& key);
    void insert(K key, V value);
    long insert_and_get_count(K key, V value);
    size_t size();
    size_t bucket_count();
    K* keys();
    V* values();
    void clear();
    void print_buckets();
};


/**
 * @class ParallelHashMap ParallelHashMap.h "src/ParallelHashMap.h"
 * @brief A thread-safe hash map supporting insertion and lookup operations.
 * @details The ParallelHashMap class is built on top of the FixedHashMap
 *    class, supporting insertion and lookup operations but not deletion as
 *    deletion is not needed in the OpenMOC application. This hash table uses
 *    chaining for collisions, as defined in FixedHashMap. It offers lock
 *    free lookups in O(1) time on average and fine-grained locking for
 *    insertions in O(1) time on average as well. Resizing is conducted
 *    periodically during inserts, although the starting table size can be
 *    chosen to limit the number of resizing operations.
 */
template <class K, class V>
class ParallelHashMap {

  /* padded pointer to hash table to avoid false sharing */
  struct paddedPointer {
    volatile long pad_L1;
    volatile long pad_L2;
    volatile long pad_L3;
    volatile long pad_L4;
    volatile long pad_L5;
    volatile long pad_L6;
    volatile long pad_L7;
    FixedHashMap<K,V> volatile* value;
    volatile long pad_R1;
    volatile long pad_R2;
    volatile long pad_R3;
    volatile long pad_R4;
    volatile long pad_R5;
    volatile long pad_R6;
    volatile long pad_R7;
    volatile long pad_R8;
  };

  private:
    FixedHashMap<K,V> *_table;
    paddedPointer *_announce;
    size_t _num_threads;
    size_t _N;
    omp_lock_t * _locks;
    size_t _num_locks;
    bool _fixed_size;
    void resize();

  public:
    ParallelHashMap(size_t M = 64, size_t L = 64);
    virtual ~ParallelHashMap();
    bool contains(K& key);
    V at(K& key);
    void update(K& key, V value);
    void insert(K key, V value);
    long insert_and_get_count(K key, V value);
    size_t size();
    size_t bucket_count();
    size_t num_locks();
    K* keys();
    V* values();
    void clear();
    void setNumThreads(int num_threads);
    void setFixedSize();
    void print_buckets();
    void realloc(size_t M);
};


/**
 * @brief Constructor initializes fixed-size table of buckets filled with empty
 *      linked lists.
 * @details The constructor initializes a fixed-size hash map with the size
 *      as an input parameter. If no size is given the default size (64)
 *      is used. Buckets are filled with empty linked lists presented as
 *      NULL pointers.
 * @param M size of fixed hash map
 */
template <class K, class V>
FixedHashMap<K,V>::FixedHashMap(size_t M) {

  /* ensure M is a power of 2 */
  if ((M & (M-1)) != 0) {
    /* if not, round up to nearest power of 2 */
    M--;
    for (size_t i = 1; i < 8 * sizeof(size_t); i*=2)
      M |= M >> i;
    M++;
  }

  /* allocate table */
  _M = M;
  _N = 0;
  _buckets = new node*[_M]();
}


/**
 * @brief Destructor deletes all nodes in the linked lists associated with each
 *      bucket in the fixed-size table and their pointers.
 */
template <class K, class V>
FixedHashMap<K,V>::~FixedHashMap() {
  /* for each bucket, scan through linked list and delete all nodes */
  for (size_t i=0; i<_M; i++) {
    node *iter_node = _buckets[i];
    while (iter_node != NULL) {
      node *next_node = iter_node->next;
      delete iter_node;
      iter_node = next_node;
    }
  }

  /* delete all buckets (now pointers to empty linked lists) */
  delete [] _buckets;
}


/**
 * @brief Determine whether the fixed-size table contains a given key.
 * @details The linked list in the bucket associated with the key is searched
 *       to determine whether the key is present.
 * @param key key to be searched
 * @return boolean value referring to whether the key is contained in the map
 */
template <class K, class V>
bool FixedHashMap<K,V>::contains(K& key) {

  /* get hash into table assuming M is a power of 2, using fast modulus */
  size_t key_hash = std::hash<K>()(key) & (_M-1);

  /* search corresponding bucket for key */
  node *iter_node = _buckets[key_hash];
  while (iter_node != NULL) {
    if (iter_node->key == key)
      return true;
    else
      iter_node = iter_node->next;
  }
  return false;
}


/**
 * @brief Determine the value associated with a given key in the fixed-size
 *      table.
 * @details The linked list in the bucket associated with the key is searched
 *      and once the key is found, the corresponding value is returned.
 *      An exception is thrown if the key is not present in the map.
 * @param key key whose corresponding value is desired
 * @return value associated with the given key
 */
template <class K, class V>
V& FixedHashMap<K,V>::at(K& key) {

  /* get hash into table assuming M is a power of 2, using fast modulus */
  size_t key_hash = std::hash<K>()(key) & (_M-1);

  /* search bucket for key and return the corresponding value if found */
  node *iter_node = _buckets[key_hash];
  while (iter_node != NULL) {
    if (iter_node->key == key)
      return iter_node->value;
    else
      iter_node = iter_node->next;
  }

  /* after the bucket has been completely searched without finding the key,
     print an error message */
  log_printf(ERROR, "Key not present in map");

  /* Should never be reached, to silence a compilation warning */
  return _buckets[0]->value;
}


/**
 * @brief Inserts a key/value pair into the fixed-size table.
 * @details The specified key value pair is inserted into the fixed-size table.
 *      If the key already exists in the table, the pair is not inserted
 *      and the function returns.
 * @param key key of the key/value pair to be inserted
 * @param value value of the key/value pair to be inserted
 */
template <class K, class V>
void FixedHashMap<K,V>::insert(K key, V value) {

  /* get hash into table using fast modulus */
  size_t key_hash = std::hash<K>()(key) & (_M-1);

  /* check to see if key already exists in map */
  if (contains(key))
    return;

  /* create new node */
  node *new_node = new node(key, value);

  /* find where to place element in linked list */
  node **iter_node = &_buckets[key_hash];
  while (*iter_node != NULL)
    iter_node = &(*iter_node)->next;

  /* place element in linked list */
  *iter_node = new_node;

  /* increment counter */
#pragma omp atomic
  _N++;
}


/**
 * @brief Inserts a key/value pair into the fixed-size table and returns the
 *      order number with which it was inserted.
 * @details The specified key value pair is inserted into the fixed-size table.
 *      If the key already exists in the table, the pair is not inserted
 *      and the function returns -1.
 * @param key key of the key/value pair to be inserted
 * @param value value of the key/value pair to be inserted
 * @return order number in which key/value pair was inserted, -1 is returned if
 *      key was already present in map.
 */
template <class K, class V>
long FixedHashMap<K,V>::insert_and_get_count(K key, V value) {

  /* get hash into table using fast modulus */
  size_t key_hash = std::hash<K>()(key) & (_M-1);

  /* check to see if key already exists in map */
  if (contains(key))
    return -1;

  /* create new node */
  node *new_node = new node(key, value);

  /* find where to place element in linked list */
  node **iter_node = &_buckets[key_hash];
  while (*iter_node != NULL)
    iter_node = &(*iter_node)->next;

  /* place element in linked list */
  *iter_node = new_node;

  /* increment counter and return number */
  size_t N;
#pragma omp atomic capture
    N = _N++;

  return (long) N;
}


/**
 * @brief Returns the number of key/value pairs in the fixed-size table.
 * @return number of key/value pairs in the map
 */
template <class K, class V>
size_t FixedHashMap<K,V>::size() {
  return _N;
}


/**
 * @brief Returns the number of buckets in the fixed-size table.
 * @return number of buckets in the map
 */
template <class K, class V>
size_t FixedHashMap<K,V>::bucket_count() {
  return _M;
}


/**
 * @brief Returns an array of the keys in the fixed-size table.
 * @details All buckets are scanned in order to form a list of all keys
 *      present in the table and then the list is returned. WARNING: The user
 *      is responsible for freeing the allocated memory once the array is no
 *      longer needed.
 * @return an array of keys in the map whose length is the number of key/value
 *      pairs in the table.
*/
template <class K, class V>
K* FixedHashMap<K,V>::keys() {

  /* allocate array of keys */
  K *key_list = new K[_N];

  /* fill array with keys */
  size_t ind = 0;
  for (size_t i=0; i<_M; i++) {
    node *iter_node = _buckets[i];
    while (iter_node != NULL) {
      key_list[ind] = iter_node->key;
      iter_node = iter_node->next;
      ind++;
    }
  }
  return key_list;
}


/**
 * @brief Returns an array of the values in the fixed-size table.
 * @details All buckets are scanned in order to form a list of all values
 *      present in the table and then the list is returned. WARNING: The user
 *      is responsible for freeing the allocated memory once the array is no
 *      longer needed.
 * @return an array of values in the map whose length is the number of
 *      key/value pairs in the table.
*/
template <class K, class V>
V* FixedHashMap<K,V>::values() {

  /* allocate array of values */
  V *values = new V[_N];

  /* fill array with values */
  size_t ind = 0;
  for (size_t i=0; i<_M; i++) {
    node *iter_node = _buckets[i];
    while (iter_node != NULL) {
      values[ind] = iter_node->value;
      iter_node = iter_node->next;
      ind++;
    }
  }
  return values;
}


/**
 * @brief Clears all key/value pairs form the hash table.
 */
template <class K, class V>
void FixedHashMap<K,V>::clear() {

  /* for each bucket, scan through linked list and delete all nodes */
  for (size_t i=0; i<_M; i++) {
    node *iter_node = _buckets[i];
    while (iter_node != NULL) {
      node *next_node = iter_node->next;
      delete iter_node;
      iter_node = next_node;
    }
  }

  /* reset each bucket to null */
  for (size_t i=0; i<_M; i++)
    _buckets[i] = NULL;

  /* reset the number of entries to zero */
  _N = 0;
}


/**
 * @brief Prints the contents of each bucket to the screen.
 * @details All buckets are scanned and the contents of the buckets are
 *      printed, which are pointers to linked lists. If the pointer is NULL
 *      suggesting that the linked list is empty, NULL is printed to the
 *      screen.
 */
template <class K, class V>
void FixedHashMap<K,V>::print_buckets() {
  log_printf(NORMAL, "Printing all buckets in the hash map...");
  for (size_t i=0; i<_M; i++) {
    if (_buckets[i] == NULL)
      log_printf(NORMAL, "Bucket %d -> NULL", i);
    else
      log_printf(NORMAL, "Bucket %d -> %p", i, _buckets[i]);
  }
}


/**
 * @brief Constructor generates initial underlying table as a fixed-sized
 *      hash map and initializes concurrency structures.
 */
template <class K, class V>
ParallelHashMap<K,V>::ParallelHashMap(size_t M, size_t L) {

  /* allocate table */
  _table = new FixedHashMap<K,V>(M);
  _fixed_size = false;

  /* get number of threads and create concurrency structures */
  _num_threads = omp_get_max_threads();
  _num_locks = L;
  _locks = new omp_lock_t[_num_locks];
  for (size_t i=0; i<_num_locks; i++)
    omp_init_lock(&_locks[i]);

  _announce = new paddedPointer[_num_threads];
  for (size_t t=0; t<_num_threads; t++)
    _announce[t].value = NULL;
}


/**
 * @brief Destructor frees memory associated with fixed-sized hash map and
 *      concurrency structures.
 */
template <class K, class V>
ParallelHashMap<K,V>::~ParallelHashMap() {
  delete _table;
  delete [] _locks;
  delete [] _announce;
}


/**
 * @brief Determine whether the parallel hash map contains a given key.
 * @details First the thread accessing the table announces its presence and
 *      which table it is reading. Then the linked list in the bucket
 *      associated with the key is searched without setting any locks
 *      to determine whether the key is present. When the thread has
 *      finished accessing the table, the announcement is reset to NULL.
 *      The announcement ensures that the data in the map is not freed
 *      during a resize until all threads have finished accessing the map.
 * @param key key to be searched
 * @return boolean value referring to whether the key is contained in the map
 */
template <class K, class V>
bool ParallelHashMap<K,V>::contains(K& key) {

  /* get thread ID */
  size_t tid = 0;
  tid = omp_get_thread_num();

  /* get pointer to table, announce it will be searched,
     and ensure consistency */
  FixedHashMap<K,V> *table_ptr;
  do {
    table_ptr = _table;
    _announce[tid].value = table_ptr;
#pragma omp flush
  } while (table_ptr != _table);

  /* see if current table contains the thread */
  bool present = table_ptr->contains(key);

  /* reset table announcement to not searching */
  _announce[tid].value = NULL;

  return present;
}


/**
 * @brief Determine the value associated with a given key.
 * @details This function follows the same algorithm as "contains" except that
 *      the value associated with the searched key is returned.
 *      First the thread accessing the table acquires the lock corresponding
 *      with the associated bucket based on the key. Then the linked list
 *      in the bucket is searched for the key. An exception is thrown if the
 *      key is not found. When the thread has finished accessing the table,
 *      it releases the lock.
 * @param key key to be searched
 * @return value associated with the key
 */
template <class K, class V>
V ParallelHashMap<K,V>::at(K& key) {

  /* If the size is fixed, simply return the value from the fixed hash map */
  if (_fixed_size)
    return _table->at(key);

  /* get thread ID */
  size_t tid = 0;
  tid = omp_get_thread_num();

  /* get pointer to table, announce it will be searched */
  FixedHashMap<K,V> *table_ptr;
  do {
    table_ptr = _table;
    _announce[tid].value = table_ptr;
#pragma omp flush
  } while (table_ptr != _table);

  /* get value associated with the key in the underlying table */
  V value = table_ptr->at(key);

  /* reset table announcement to not searching */
  _announce[tid].value = NULL;

  return value;
}


/**
 * @brief Insert a given key/value pair into the parallel hash map.
 * @details First the underlying table is checked to determine if a resize
 *      should be conducted. Then, the table is checked to see if it
 *      already contains the key. If so, the key/value pair is not inserted
 *      and the function returns. Otherwise, the lock of the associated
 *      bucket is acquired and the key/value pair is added to the bucket.
 * @param key key of the key/value pair to be inserted
 * @param value value of the key/value pair to be inserted
 */
template <class K, class V>
void ParallelHashMap<K,V>::insert(K key, V value) {
  /* check if resize needed */
  if (2*_table->size() > _table->bucket_count())
    resize();

  /* check to see if key is already contained in the table */
  if (contains(key))
    return;

  /* get lock hash */
  size_t lock_hash = (std::hash<K>()(key) & (_table->bucket_count() - 1))
    % _num_locks;

  /* acquire lock */
  omp_set_lock(&_locks[lock_hash]);

  /* insert value */
  _table->insert(key, value);

  /* release lock */
  omp_unset_lock(&_locks[lock_hash]);
}


/**
 * @brief Updates the value associated with a key in the parallel hash map.
 * @details The thread first acquires the lock for the bucket associated with
 *      the key is acquired, then the linked list in the bucket is searched
 *      for the key. If the key is not found, an exception is returned. When
 *      the key is found, the value is updated and the lock is released.
 * @param key the key of the key/value pair to be updated
 * @param value the new value for the key/value pair
 */
template <class K, class V>
void ParallelHashMap<K,V>::update(K& key, V value) {

  /* get lock hash */
  size_t lock_hash = (std::hash<K>()(key) & (_table->bucket_count() - 1))
    % _num_locks;

  /* acquire lock */
  omp_set_lock(&_locks[lock_hash]);

  /* insert value */
  _table->at(key) = value;

  /* release lock */
  omp_unset_lock(&_locks[lock_hash]);
}


/**
 * @brief Insert a given key/value pair into the parallel hash map and return
      the order number.
 * @details First the underlying table is checked to determine if a resize
 *      should be conducted. Then, the table is checked to see if it
 *      already contains the key. If so, the key/value pair is not inserted
 *      and the function returns. Otherwise, the lock of the associated
 *      bucket is acquired and the key/value pair is added to the bucket.
 * @param key key of the key/value pair to be inserted
 * @param value value of the key/value pair to be inserted
 * @return order number in which the key/value pair was inserted, -1 if it
 *      already exists
 */
template <class K, class V>
long ParallelHashMap<K,V>::insert_and_get_count(K key, V value) {

  /* check if resize needed */
  if (2*_table->size() > _table->bucket_count())
    resize();

  /* check to see if key is already contained in the table */
  if (contains(key))
    return -1;

  /* get lock hash */
  size_t lock_hash = (std::hash<K>()(key) & (_table->bucket_count() - 1))
    % _num_locks;

  /* acquire lock */
  omp_set_lock(&_locks[lock_hash]);

  /* insert value */
  long N =_table->insert_and_get_count(key, value);

  /* release lock */
  omp_unset_lock(&_locks[lock_hash]);

  return N;
}


/**
 * @brief Resizes the underlying table to twice its current capacity.
 * @details In a thread-safe manner, this procedure resizes the underlying
 *    FixedHashMap table to twice its current capacity using locks and the
 *    announce array. First, all locks are set in order to block inserts and
 *    prevent deadlock. A new table is allocated of twice the size and all
 *    key/value pairs from the old table, then the pointer is switched to the
 *    new table and locks are released. Finally the memory needs to be freed.
 *    To prevent threads currently reading the table from encountering
 *    segmentation faults, the resizing threads waits for the announce array
 *    to be free of references to the old table before freeing the memory.
 */
template <class K, class V>
void ParallelHashMap<K,V>::resize() {

  /* do not resize if fixed size */
  if (_fixed_size)
    return;

  /* acquire all locks in order */
  for (size_t i=0; i<_num_locks; i++)
    omp_set_lock(&_locks[i]);

  /* recheck if resize needed */
  if (2*_table->size() <= _table->bucket_count()) {
    /* release locks */
    for (size_t i=0; i<_num_locks; i++)
      omp_unset_lock(&_locks[i]);

    return;
  }

  /* allocate new hash map of double the size */
  FixedHashMap<K,V> *new_map =
    new FixedHashMap<K,V>(2*_table->bucket_count());

  /* get keys, values, and number of elements */
  K *key_list = _table->keys();
  V *value_list = _table->values();

  /* insert key/value pairs into new hash map */
  for (size_t i=0; i<_table->size(); i++)
    new_map->insert(key_list[i], value_list[i]);

  /* save pointer of old table */
  FixedHashMap<K,V> *old_table = _table;

  /* reassign pointer */
  _table = new_map;
#pragma omp flush

  /* release all locks */
  for (size_t i=0; i<_num_locks; i++)
    omp_unset_lock(&_locks[i]);

  /* delete key and value list */
  delete [] key_list;
  delete [] value_list;

  /* wait for all threads to stop reading from the old table */
  for (size_t i=0; i<_num_threads; i++) {
    while (_announce[i].value == old_table) {
#pragma omp flush
      continue;
    }
  }

  /* free memory associated with old table */
  delete old_table;
}


/**
 * @brief Returns the number of key/value pairs in the underlying table.
 * @return number of key/value pairs in the map
 */
template <class K, class V>
size_t ParallelHashMap<K,V>::size() {
  return _table->size();
}


/**
 * @brief Returns the number of buckets in the underlying table.
 * @return number of buckets in the map
 */
template <class K, class V>
size_t ParallelHashMap<K,V>::bucket_count() {
  return _table->bucket_count();
}


/**
 * @brief Returns the number of locks in the parallel hash map.
 * @return number of locks in the map
 */
template <class K, class V>
size_t ParallelHashMap<K,V>::num_locks() {
  return _num_locks;
}


/**
 * @brief Returns an array of the keys in the underlying table.
 * @details All buckets are scanned in order to form a list of all keys
 *      present in the table and then the list is returned. Threads
 *      announce their presence to ensure table memory is not freed
 *      during access. WARNING: The user is responsible for freeing the
 *      allocated memory once the array is no longer needed.
 * @return an array of keys in the map whose length is the number of key/value
 *      pairs in the table.
 */
template <class K, class V>
K* ParallelHashMap<K,V>::keys() {

  /* get thread ID */
  size_t tid = 0;
  tid = omp_get_thread_num();

  /* get pointer to table, announce it will be searched */
  FixedHashMap<K,V> *table_ptr;
  do {
    table_ptr = _table;
    _announce[tid].value = table_ptr;
#pragma omp flush
  } while (table_ptr != _table);

  /* get key list */
  K* key_list = table_ptr->keys();

  /* reset table announcement to not searching */
  _announce[tid].value = NULL;

  return key_list;
}


/**
 * @brief Returns an array of the values in the underlying table.
 * @details All buckets are scanned in order to form a list of all values
 *      present in the table and then the list is returned. Threads
 *      announce their presence to ensure table memory is not freed
 *      during access. WARNING: The user is responsible for freeing the
 *      allocated memory once the array is no longer needed.
 * @return an array of values in the map whose length is the number of key/value
 *      pairs in the table.
 */
template <class K, class V>
V* ParallelHashMap<K,V>::values() {

  /* get thread ID */
  size_t tid = 0;
  tid = omp_get_thread_num();

  /* get pointer to table, announce it will be searched */
  FixedHashMap<K,V> *table_ptr;
  do {
    table_ptr = _table;
    _announce[tid].value = table_ptr;
#pragma omp flush
  } while (table_ptr != _table);

  /* get value list */
  V *value_list = table_ptr->values();

  /* reset table announcement to not searching */
  _announce[tid].value = NULL;

  return value_list;
}


/**
 * @brief Re-allocate the threaded _announce vector to the desired number of
 *        threads.
 * @param num_threads the desired number of threads
 */
template <class K, class V>
void ParallelHashMap<K,V>::setNumThreads(int num_threads) {

  _num_threads = num_threads;

  /* Allocate for the desired number of threads */
  delete[] _announce;
  _announce = new paddedPointer[_num_threads];
  for (size_t t=0; t<_num_threads; t++)
    _announce[t].value = NULL;
}


/**
 * @brief Prevents the parallel hash map from further resizing.
 */
template <class K, class V>
void ParallelHashMap<K,V>::setFixedSize() {
  _fixed_size = true;
}


/**
 * @brief Clears all key/value pairs form the hash table.
 */
template <class K, class V>
void ParallelHashMap<K,V>::clear() {

  /* acquire all locks in order */
  for (size_t i=0; i<_num_locks; i++)
    omp_set_lock(&_locks[i]);

  /* clear underlying fixed table */
  _table->clear();

  /* release all locks in order */
  for (size_t i=0; i<_num_locks; i++)
    omp_unset_lock(&_locks[i]);
}


/**
 * @brief Prints the contents of each bucket to the screen.
 * @details All buckets are scanned and the contents of the buckets are
 *      printed, which are pointers to linked lists. If the pointer is NULL
 *      suggesting that the linked list is empty, NULL is printed to the
 *      screen. Threads announce their presence to ensure table memory is
 *      not freed during access.
 */
template <class K, class V>
void ParallelHashMap<K,V>::print_buckets() {

  /* get thread ID */
  size_t tid = 0;
  tid = omp_get_thread_num();

  /* get pointer to table, announce it will be searched */
  FixedHashMap<K,V> *table_ptr;
  do {
    table_ptr = _table;
    _announce[tid].value = table_ptr;
#pragma omp flush
  } while (table_ptr != _table);

  /* print buckets */
  table_ptr->print_buckets();

  /* reset table announcement to not searching */
  _announce[tid].value = NULL;
}

/**
 * @brief Reallocates the underlying table to the desired size.
 * @param M The requested new size
 */
template <class K, class V>
void ParallelHashMap<K,V>::realloc(size_t M) {

  /* delete old table */
  delete _table;

  /* allocate new table */
  _table = new FixedHashMap<K,V>(M);
}
#endif
