/*
	A header-based library for Mantissa (Significand) segmentation.
	Author: harunadess

	Based on the basic idea of mantissa segmentation as detailed in "A Customized Precision format based on Mantissa Segmentation" (https://doi.org/10.1002/cpe.5418)

	Provides classes for the adoption of mantissa segmentation as a techinque into applications, specifically those
	that are memory bound (and especially iterative algorithms).
	The basic concept is that an algorithm may perform some initial calculations in reduced precision (i.e. using only some of the mantissa)
	while still having access to the full exponent of a double precision number, and without the need for data duplication in multiple precisions.
	At some later specified point, the switch is then made to full precision to obtain a more accurate solution.
	This can be done in-place by reading more segments, or by copying data to a double precision array (not in-place) but achieving better
	full precision performance.
	
	Copyright (c) 2020 harunadess

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
*/

#ifndef __MANSEG_LIB_H__
#define __MANSEG_LIB_H__

#include <stdint.h>
#include <immintrin.h>

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef uint_fast64_t doublerep;

    /* Highest achievable precision with a single segment - i.e. TwoSegArray<false> */
    constexpr double MaxSingleSegmentPrecision = 1e-6;
    /* decimal precision: num_mantissa_bits*log10(2) = 6~ */
	/* so we have approximately 6 digits of decimal precision */
    constexpr double AdaptivePrecisionBound = 5e-5;

    /*
        Class representing the "head" segment of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<false>.
        It contains the sign bit, full 11 bit exponent, and 20 bits of mantissa, for a total
        of 32 bits.
        This gives less precision than IEEE-754 standard float.
        It has a maximum precision of roughly 1e-6, with a recommended precision bound of 1.5e-6.
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
            return head_vd[0];
        }

        float* head;
    };

    /*
        Class representing the "pair" of segments of a double.
        It is defined for ease of manipulating the values in an object of type TwoSegArray<true>.
        Each segment (which consists of a head and tail) is 32 bits.
        The "head" contains the sign bit, full 11-bit exponent, and 20 bits of mantissa, 
        for a total of 32 bits.
        The "tail" contains the remaining 32 bits of mantissa.
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

		Pair(const Head& o)
		{
			*head = *o.head;
			*tail = 0;
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
            return seg_vd[0];
        }

        float* head;
        float* tail;
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

		Note: the class deconstructor does not free space automatically, so use the delSegments and del functions to clean up.
    */
    template<>
    class TwoSegArray<true>
    {
    public:
        TwoSegArray() 
		{
			heads = nullptr;
			tails = nullptr;
		}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new float[length];
            tails = new float[length] (); // initially zero tails array
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray() { }

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
            Returns object of the same type, with pointers to the same data values
			as this object.
        */
        TwoSegArray<true> createFullPrecision()
        {
            return TwoSegArray<true>(heads, tails);
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new float[length];
            tails = new float[length] ();
        }

		bool isAlloc() { return (heads != nullptr) && (tails != nullptr); }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del()
        {
            if(heads != nullptr) delete[] heads;
            if(tails != nullptr) delete[] tails;
			heads = nullptr;
			tails = nullptr;
        }

    private:
        float* heads;
        float* tails;
    };

    /*
        Specialisation of TwoSegArray.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
		Note: the class deconstructor does not free space automatically, so use the delSegments and del functions to clean up.
    */
    template<>
    class TwoSegArray<false>
    {
    public:
        TwoSegArray()
		{
			heads = nullptr;
			tails = nullptr;
		}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new float[length];
            tails = new float[length] (); // initially zero tails array
        }

        TwoSegArray(float* heads, float* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray() { }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
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
			tails[id] = seg_v[0];
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
            heads = new float[length];
            // we should zero tails when allocating heads array for
            // to avoid unexpected behaviour
            tails = new float[length] ();
        }

		bool isAlloc() { return (heads != nullptr) && (tails != nullptr); }

        /*
            Deletes the values of the dynamic arrays used to store values in
            the array.
            NOTE: this should only be called by one object with references to
            the same set of values (such as object created using createFullPrecision).
        */
        void del()
        {
            if(heads != nullptr) delete[] heads;
            if(tails != nullptr) delete[] tails;   
			heads = nullptr;
			tails = nullptr;
        }

    private:
        float* heads;
        float* tails;
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
        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T&& other) noexcept
    {
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
		convenience typename for TwoSegArray that only accesses the "head" segment 
		of a double precision number
		i.e. reduced precision
	*/
    using HeadsArray = TwoSegArray<false>;

	/* 
		convenience typename for TwoSegArray that accesses the full 64 bits of
		a double precision number
		i.e. double precision
	*/
    using PairsArray = TwoSegArray<true>;

    /*
        Convenient type for use of TwoSegArray<false> and TwoSegArray<true> without having to manage two separate sets of arrays.
		Note: the class deconstructor does not free space automatically, so use the delSegments and del functions to clean up.
        
        @heads - used to access only the upper 32 bits of a double : [sign(1), exp(11), mantissa(20)]
        @pairs - used to access all 64 bits of a double : [sign(1), exp(11), mantissa(52)] (slow, but does not require extra memory)
        @full - used to access all 64 bits of double, in the standard IEEE method (fast, recommended, but requires extra space)
    */
    class ManSegArray
    {
    public:
        HeadsArray heads; 		// provides access to values in reduced precision (i.e. the upper 32 bits of a double precision number)
        PairsArray pairs; 		// provides access to values in full precision, by converting values from segments (i.e. all 64 bits of a double precision number)
        double* full; 			// standard double precision access; requires use of copyToIEEEdouble, or manual population
        uint_fast64_t length; 	// length of allocated array; should be manually set if parameterised constructor/alloc is not used.

        ManSegArray() { full = nullptr; length = 0; }

        ~ManSegArray() { }

        ManSegArray(const uint_fast64_t& length)
        {
            this->length = length;
            heads.alloc(length);
            pairs = heads.createFullPrecision();
			full = nullptr;
        }

        /*
            Allocates length elements to the array dynamically.
            Note: this array space is used for both heads and pairs
        */
        void alloc(const uint_fast64_t& length)
        {
			this->length = length;
            heads.alloc(length);
            pairs = heads.createFullPrecision();
        }

        /*
            Implements precision switching by allocating length space for copying the full 64-bit values from
			the pairs array to the full doubles array.

            This is implemented using omp parallel, however it can also be accomplished by the user, as full is
			publically available. Note: the length variable should be set if a user implemented copy is performed.
        */
        void copytoIEEEdouble()
        {
            full = new double[length];

            #pragma omp parallel for
			for(int i = 0; i < length; ++i)
				full[i] = heads[i];
        }

        /* 
            Deletes space allocated to the segments arrays.
            WARNING: should only be called once, as heads and pairs share the array space.
        */
        void delSegments() { heads.del(); if(full == nullptr) length = 0; }

		/*
			Deletes space allocated to full IEEE double precision array.
			Must be called in order to free space.
		*/
		void del() { if(full != nullptr) delete[] full; full = nullptr; if(!heads.isAlloc()) length = 0; }
    };
}

#endif
