// not sure if this right
// I've never written a .h file before

#if HAVEFP
double drand48();
#else
long irand48(m);
long krand48(xsubi, m);
#endif

long lrand48();
long mrand48();
static void next();

void srand48(seedval);

unsigned short * seed48(seed16v);

void lcong48(param);
