#pragma once
#ifdef _UNIX
#include <pthread.h>
namespace { struct X { X() { pthread_getconcurrency(); } } x; }
#endif
