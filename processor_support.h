#ifndef PROCESSOR_SUPPORT_H_
#define PROCESSOR_SUPPORT_H_

// Utility class ProcessorSupport provides POPCNTenabled() to determine
// processor support for POPCNT instruction. It uses CPUID to
// retrieve the processor capabilities.
// for Intel ICC compiler __cpuid() is an intrinsic 
// for Microsoft compiler __cpuid() is provided by #include <intrin.h>
// for GCC compiler __get_cpuid() is provided by #include <cpuid.h>

// Intel compiler defines __GNUC__, so this is needed to disambiguate

// CPUID / POPCNT capability detection is x86-only. On other architectures
// (e.g. ARM aarch64) there is no cpuid.h and no __cpuid(); POPCNTenabled()
// returns false and callers use the portable popcount fallback.
#if defined(__amd64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#   define PROCESSOR_IS_X86
#endif

#if defined(PROCESSOR_IS_X86)
#   if defined(__INTEL_COMPILER)
#       define USING_INTEL_COMPILER
#   elif defined(__GNUC__)
#       define USING_GCC_COMPILER
#       include <cpuid.h>
#   elif defined(_MSC_VER)
// __MSC_VER defined by Microsoft compiler
#       define USING_MSC_COMPILER
#   endif
#endif

struct regs_t {unsigned int EAX, EBX, ECX, EDX;};
#define BIT(n) ((1<<n))

class ProcessorSupport {

#ifdef POPCNT_CAPABILITY 

public: 
    ProcessorSupport() { } 
    bool POPCNTenabled()
    {
    // from: Intel® 64 and IA-32 Architectures Software Developer’s Manual, 325462-036US,March 2013
    //Before an application attempts to use the POPCNT instruction, it must check that the
    //processor supports SSE4.2
    //“(if CPUID.01H:ECX.SSE4_2[bit 20] = 1) and POPCNT (if CPUID.01H:ECX.POPCNT[bit 23] = 1)”
    //
    // see p.272 of http://download.intel.com/products/processor/manual/253667.pdf available at
    // http://www.intel.com/content/www/us/en/processors/architectures-software-developer-manuals.html
    // Also http://en.wikipedia.org/wiki/SSE4 talks about available on Intel & AMD processors

    regs_t regs;

    try {
#if ( defined(USING_INTEL_COMPILER) || defined(USING_MSC_COMPILER) )
        __cpuid((void *) &regs,0); // test if __cpuid() works, if not catch the exception
        __cpuid((void *) &regs,0x1); // POPCNT bit is bit 23 in ECX
#elif defined(USING_GCC_COMPILER)
        __get_cpuid(0x1, &regs.EAX, &regs.EBX, &regs.ECX, &regs.EDX);
#else
        // No CPUID available (e.g. non-x86 architecture such as ARM): report
        // POPCNT as unsupported so the portable popcount fallback is used.
        return false;
#endif
        if( !( (regs.ECX & BIT(20)) && (regs.ECX & BIT(23)) ) ) return false;
    }
    catch (int e) {
        return false;
    }
    return true;
    }

#endif // POPCNT_CAPABILITY
};

#endif /*PROCESSOR_SUPPORT_H_*/




