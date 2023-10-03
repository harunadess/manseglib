/*
	A header-based library for Mantissa (Significand) segmentation.
	Author: harunadess

	Based on the basic idea of mantissa segmentation as detailed in "A Customized Precision format based on Mantissa Segmentation" (https://doi.org/10.1002/cpe.5418)
	
	Debug/trial version of the float/vector version of the library.
*/

#ifndef __MANSEG_LIB_H__
#define __MANSEG_LIB_H__

#include <stdint.h>
#include <immintrin.h>
#include <iostream>

#ifndef DEBUG_FUNC_MANSEG
#define DEBUG_FUNC_MANSEG
// print binary of double
void printBinary(double d)
{
    unsigned long l = *reinterpret_cast<unsigned long*>(&d);
    for(int i = 63; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 63 || i == 52)
            std::cout << " ";
        else if(i == 32)
            std::cout << "|";
    }
    std::cout << std::endl;
}

// print binary of float
void printBinary(float f)
{
    unsigned int l = *reinterpret_cast<unsigned int*>(&f);
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i == 31 || i == 23)
            std::cout << " ";
    }
    std::cout << std::endl;
}

// print binary of int
void printBinary(int l)
{
    for(int i = 31; i >= 0; --i)
    {
        std::cout << ((l >> i) & 1);
        if(i % 8 == 0 && i > 0)
            std::cout << " ";
    }
}
#endif

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef uint_fast64_t doublerep;

    /* Highest achievable precision with a single segment - i.e. TwoSegArray<false> */
    constexpr double MaxSingleSegmentPrecision = 1e-5;
    /* decimal precision: num_mantissa_bits*log10(2) = 6.02... */
    constexpr double AdaptivePrecisionBound = 5e-5;

    /*
        Class representing the "head" segment of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<false>.
        It contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, for a total
        of 32 bits.
        This gives less precision than IEEE-754 standard float.
        It has a maximum precision of roughly 1e-5, with a recommended precision bound of 5e-5.
    */
    class Head
    {
    public:
        Head(float* head)
            :head(head)
        {}

        Head(const Head& o)
        {
            *head = *o.head;
        }

        /* 
            These equals operators must either be combined into a generic template function and
            therefore could result in extra conversions when copying data.
            If we want to keep the specialisation, then we need to either move the definition into
            a separate .cpp file and compile a .o file for the library, or we can declare them
            inline, as that allows multiple redefinitions.
            Considering it is only for a few functions, (the TwoSegArray class specialisations must be
            kept in this file) it seems like a reasonable compromise.
        */
        template<typename T>
        inline Head& operator=(const T& rhs);
        template<typename T>
        inline Head& operator=(const T&& other) noexcept;
        
        template<typename T>
        double operator+=(const T& rhs);
        template<typename T>
        double operator-=(const T& rhs);
        template<typename T>
        double operator*=(const T& rhs);
        template<typename T>
        double operator/=(const T& rhs);

        operator double() const
        {
            __m128 head_v = _mm_set_ps(0.0f, 0.0f, *head, 0.0f);
            __m128d head_vd = _mm_castps_pd(head_v);
            double d;
            _mm_store_sd(&d, head_vd);
            return d;
        }

        float* head;
        static constexpr int segmentBits = 32;
        static constexpr __m128i head_mask = {(int_fast64_t)(0xFFFFFFFF00000000), (int_fast64_t)(0x0000000000000000)};
    };

    /*
        Class representing the "pair" of segments of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<true>.
        Each segment (which consists of a head and tail) is 32 bits.
        The head contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, 
        for a total of 32 bits.
        The tail segment contains the remaining 32 bits of mantissa.
    */
    class Pair
    {
    public:
        Pair(float* head, float* tail)
            :head(head), tail(tail)
        {}

        Pair(const Pair& o)
        {
            *head = *o.head;
            *tail = *o.tail;
        }

        template<typename T>
        inline Pair& operator=(const T& rhs);
        template<typename T>
        inline Pair& operator=(const T&& other) noexcept;

        template<typename T>
        double operator+=(const T& rhs);
        template<typename T>
        double operator-=(const T& rhs);
        template<typename T>
        double operator*=(const T& rhs);
        template<typename T>
        double operator/=(const T& rhs);

        operator double() const
        {
            __m128 seg_v = _mm_set_ps(0.0f, 0.0f, *head, *tail);
            __m128d seg_vd = _mm_castps_pd(seg_v);
            double d;
            _mm_store_sd(&d, seg_vd);
            return d;
        }

        float* head;
        float* tail;
        static constexpr int segmentBits = 32;
        static constexpr __m128i head_mask = {(int_fast64_t)(0xFFFFFFFF00000000), (int_fast64_t)(0x0000000000000000)};
        // static constexpr __m128i lower_64_mask = {(int_fast64_t)(0xFFFFFFFFFFFFFFFF), (int_fast64_t)(0x0000000000000000)};
    };

    /*
        Array type class for performing operations on double values, conceptually split into two 32 bit segments - "head" and "tail".
        The user is required to manually clean up memory after allocation using the del() function, as the deconstructor does not handle this. This is due to the
        way that shifting from low to high precision is handled, by just reading the extra segments in the <true> variant.
        
        <false> specialisation - Operations performed on values in the array only modify the "head" segment (i.e. the first 32 bits of a double value) unless specified otherwise.
        <true> specialisation - Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<bool useTail>
    class TwoSegArray; // explicitly specialise this below.

    /*
        Specialisation of TwoSegArray.
        User is required to manage de-allocation of memory manually, using the del()
        function, as the deconstructor does not free this space.
        Operations performed on values in the array are exactly the same as standard IEEE-754 doubles, but a bit slower due to combining of head and tail segments.
    */
    template<>
    class TwoSegArray<true>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const uint_fast64_t& length)
        {
            // allocate one element extra
            heads = new float[1 + length];
            tails = new float[1 + length] (); // initially zero tails array
            // following two lines put the extra element of the array at the front
            // i.e. position 0
            // this makes vector operations easier to deal with.
            ++heads;
            ++tails;
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        Pair operator[](const uint_fast64_t& id)
        {
            return Pair(&heads[id], &tails[id]);
        }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            __m128d d_v = _mm_set_pd(0.0, d);
            __m128 seg_v = _mm_castpd_ps(d_v);
            tails[id] = seg_v[0];
            heads[id] = seg_v[1];
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            __m128d d_v = _mm_set_pd(0.0, d);
            __m128 seg_v = _mm_castpd_ps(d_v);
            tails[id] = seg_v[0];
            heads[id] = seg_v[1];
        }

        double read(const uint_fast64_t& id)
        {
            return static_cast<double>(Pair(&heads[id], &tails[id]));
        }

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegArray<true> createFullPrecision()
        {
            return *this;
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new float[1 + length];
            tails = new float[1 + length] ();

            ++heads;
            ++tails;
        }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */

        void del()
        {
            if(heads != nullptr) delete[] (heads - 1);
            if(tails != nullptr) delete[] (tails - 1);
        }

    private:
        float* heads;
        float* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };

    /*
        Specialisation of TwoSegArray.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
    */
    template<>
    class TwoSegArray<false>
    {
    public:
        TwoSegArray() {}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new float[1 + length];
            tails = new float[1 + length] (); // initially zero tails array

            ++heads;
            ++tails;
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray()
        {
            heads = nullptr;
            tails = nullptr;
        }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            // double d = t;
            // float *pd = (float*)&d;
            // __m128 d_v = _mm_load_ss(pd + 1);
            // _mm_store_ss(&heads[id], d_v);
            double d = t;
            __m128d d_v = _mm_set_pd(0.0, d);
            __m128 seg_v = _mm_castpd_ps(d_v);
            heads[id] = seg_v[1];
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            __m128d d_v = _mm_set_pd(0.0, d);
            __m128 seg_v = _mm_castpd_ps(d_v);
            heads[id] = seg_v[1];
        }

        double read(const uint_fast64_t& id)
        {
            return static_cast<double>(Head(&heads[id]));
        }

        Head operator[](const uint_fast64_t& id)
        {
            return Head(&heads[id]);
        }

        /*
            Increases the precision of the operations performed by returning an
            object of type TwoSegArray<true> that has pointers to the values
            in the existing array.
        */
        TwoSegArray<true> createFullPrecision()
        {
            return TwoSegArray<true>(heads, tails);
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new float[1 + length];
            // we should zero tails when allocating heads array for
            // to avoid unexpected behaviour
            tails = new float[1 + length] ();

            ++heads;
            ++tails;
        }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del()
        {
            if(heads != nullptr) delete[] (heads - 1);
            if(tails != nullptr) delete[] (tails - 1);   
        }

    private:
        float* heads;
        float* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;
    };


    /**
     * Head functions
    */

    template<>
    inline Head& Head::operator=(const Head& other)
    {
        *head = *other.head;
        return *this;
    }

    template<>
    inline Head& Head::operator=(const Head&& other) noexcept
    {
        *head = *other.head;
        return *this;
    }

    template<>
    inline Head& Head::operator=(const Pair& other)
    {
        *head = *other.head;
        return *this;
    }

    template<>
    inline Head& Head::operator=(const Pair&& other) noexcept
    {
        *head = *other.head;
        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T& other)
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        *head = seg_v[1];

        // double d = other;
        // __m128 f_v = _mm_loadu_ps(reinterpret_cast<float*>(&d) + 1);
        // _mm_store_ss(&head, f_v);

        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T&& other) noexcept
    {
        // double d = other;
        // float *pd = (float*)&d;
        // __m128 d_v = _mm_load_ss(pd + 1);
        // _mm_store_ss(head, d_v);
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        *head = seg_v[1];

        return *this;
    }

    template<typename T>
    double Head::operator+=(const T& rhs)
    {
        double t = *this;
        t += rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator+(Head lhs, const T& rhs)
    {
        double d = lhs;
        d += rhs;
        return d;
    }

    template<typename T>
    double Head::operator-=(const T& rhs)
    {
        double t = *this;
        t -= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Head lhs, const T& rhs)
    {
        double d = lhs;
        d -= rhs;
        return d;
    }

    template<typename T>
    double Head::operator*=(const T& rhs)
    {
        double t = *this;
        t *= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Head lhs, const T& rhs)
    {
        double d = lhs;
        d *= rhs;
        return d;
    }

    template<typename T>
    double Head::operator/=(const T& rhs)
    {
        double t = *this;
        t /= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Head lhs, const T& rhs)
    {
        double d = lhs;
        d /= rhs;
        return d;
    }

    /**
     * ManSegPair functions
    */
    template<>
    inline Pair& Pair::operator=(const Head& other)
    {
        *head = *other.head;
        return *this;
    }

    template<>
    inline Pair& Pair::operator=(const Head&& other) noexcept
    {
        *head = *other.head;
        return *this;
    }

    template<>
    inline Pair& Pair::operator=(const Pair& other)
    {
        *head = *other.head;
        *tail = *other.tail;
        return *this;
    }

    template<>
    inline Pair& Pair::operator=(const Pair&& other) noexcept
    {
        *head = *other.head;
        *tail = *other.tail;
        return *this;
    }

    template<typename T>
    inline Pair& Pair::operator=(const T& other)
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        *tail = seg_v[0];
        *head = seg_v[1];

        return *this;
    }

    template<typename T>
    inline Pair& Pair::operator=(const T&& other) noexcept
    {
        double d = other;
        __m128d d_v = _mm_set_pd(0.0, d);
        __m128 seg_v = _mm_castpd_ps(d_v);
        *tail = seg_v[0];
        *head = seg_v[1];

        return *this;
    }

    template <typename T>
    double Pair::operator+=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t += r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator+(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d += r;
        return d;
    }

    template<typename T>
    double Pair::operator-=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t -= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d -= r;
        return d;
    }

    template<typename T>
    double Pair::operator*=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t *= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d *= r;
        return d;
    }

    template<typename T>
    double Pair::operator/=(const T& rhs)
    {
        double t = *this;
        double r = rhs;
        t /= r;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Pair lhs, const T& rhs)
    {
        double d = lhs;
        double r = rhs;
        d /= r;
        return d;
    }

    /*
        Convenient type for use of TwoSegArray<false> and TwoSegArray<true> without having to manage two separate sets of arrays.
        
        @heads - used to access only the first 32 bits of a double : [sign(1), exp(11), mantissa(20)]
        @pairs - used to access all 64 bits of a double : [sign(1), exp(11), mantissa(52)]
    */
    class ManSegArray
    {
    public:
        TwoSegArray<false> heads;
        TwoSegArray<true> pairs;

        ManSegArray() {}

        ManSegArray(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /*
            Allocates length elements to the array dynamically.
            Note: this array space is used for both heads and pairs
        */
        inline void alloc(const uint_fast64_t& length)
        {
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /* 
            Deletes space allocated to the array.
            WARNING: should only be called once, as heads and pairs share the array space.
        */
        inline void del() { heads.del(); }
    };

    using HeadsArray = TwoSegArray<false>;
    using PairsArray = TwoSegArray<true>;
}

#endif
