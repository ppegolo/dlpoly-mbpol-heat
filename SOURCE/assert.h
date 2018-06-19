#ifndef DISABLE_ASSERT
# define assert(x) if (.not.(x)) call afailed(__FILE__,__LINE__)
#else
# define assert(x)
#endif
