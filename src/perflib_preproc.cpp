#define COMMENT(x)

#if defined(WITHPERF)
#define PERFINIT call perfinit
#define PERFON(str) call perfon(str)
#define PERFOFF call perfoff
#define C_PERFON(str,strlen) perfon(str,strlen)
#define C_PERFOFF() perfoff()
#define PERFOUT(str) call perfout(str)
#define PERF_GET(str,x,y) call perf_get(str,x,y)
#define PERF_RESET(str) call perf_reset(str)
#define PERF_CONTEXT_START(str) call perf_context_start(str)
#define PERF_CONTEXT_END call perf_context_end
#if (WITHPERF==2)
COMMENT( /*all routines are analyzed*/)
#define PERFON_I(str) call perfon(str)
#define PERFOFF_I call perfoff
#else
COMMENT(/* inner routines are not monitored to keep impact on runtime low */)
#define PERFON_I(str)
#define PERFOFF_I
#endif
#else
#define PERFINIT
#define PERFON(str)
#define PERFOFF
#define C_PERFON(str,strlen)
#define C_PERFOFF()
#define PERFON_I(str)
#define PERFOFF_I
#define PERFOUT(str)
#define PERF_GET(str,x,y)
#define PERF_RESET(str)
#define PERF_CONTEXT_START(str)
#define PERF_CONTEXT_END
#endif

#if defined(WITH_LIKWID)
#define LIKWID_INIT call likwid_markerInit()
#define LIKWID_ON(x) call likwid_markerStart(x)
#define LIKWID_OFF(x) call likwid_markerStop(x)
#define LIKWID_CLOSE call likwid_markerClose()
#else
#define LIKWID_HEADER
#define LIKWID_INIT
#define LIKWID_ON(x)
#define LIKWID_OFF(x)
#define LIKWID_CLOSE
#endif
