/*
	A header-based library for Mantissa (Significand) segmentation.
	Author: harunadess

	Based on the basic idea of mantissa segmentation as detailed in "A Customized Precision format based on Mantissa Segmentation" (https://doi.org/10.1002/cpe.5418)
	
	This is the initial attempt at creating the library, utilising integers and shifting to accomplish data conversion from double precision to segments.
	It was transformed into the vector code that was used in benchmarks, and is found in manseglib.hpp

	It is not compatible with code that uses manseglib.hpp in it's current form.
	
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

namespace ManSeg
{
    // internal representation of double, so it can be easily manipulated
    typedef uint_fast64_t doublerep;

    /* Highest achievable precision with a single segment - i.e. TwoSegArray<false> */
    constexpr double MaxSingleSegmentPrecision = 1e-5;
    /* decimal precision: num_mantissa_bits*log10(2) = 6.02... */
    constexpr double AdaptivePrecisionBound = 5e-5;

	union doublepun { double d; uint64_t l; };
	
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
        Head(uint32_t* head)
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
            doublepun conv;
			conv.l = *head;
			conv.l <<= segmentBits;
			return conv.d;
        }

        uint32_t* head;
        static constexpr int segmentBits = 32;
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
        Pair(uint32_t* head, uint32_t* tail)
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
			doublepun conv;
			conv.l = *head;
			conv.l <<= segmentBits;
			conv.l |= *tail;
			return conv.d;
        }

        uint32_t* head;
        uint32_t* tail;
        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;

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
        TwoSegArray() 
		{
			heads = nullptr;
			tails = nullptr;
		}

        TwoSegArray(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            tails = new uint32_t[length] ();
        }

        TwoSegArray(uint32_t* heads, uint32_t* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray() { }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            double d = t;
			doublepun conv;
			conv.d = d;
            tails[id] = (conv.l & tailMask);
            heads[id] = (conv.l >> segmentBits);
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
			doublepun conv;
			conv.d = d;
            tails[id] = (conv.l & tailMask);
            heads[id] = (conv.l >> segmentBits);
        }

        double read(const uint_fast64_t& id)
        {
            return static_cast<double>(Pair(&heads[id], &tails[id]));
        }

        Pair operator[](const uint_fast64_t& id)
        {
            return Pair(&heads[id], &tails[id]);
        }

        /*
            Returns *this (as we do not have any precision increase to do)
        */
        TwoSegArray<true> createFullPrecision()
        {
            return TwoSegArray<true>(heads, tails);
        }

        void alloc(const uint_fast64_t& length)
        {
            heads = new uint32_t[length];
            tails = new uint32_t[length] ();
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
        uint32_t* heads;
        uint32_t* tails;

        static constexpr int segmentBits = 32;
        static constexpr uint32_t tailMask = ~0;
        static constexpr doublerep headMask = static_cast<doublerep>(tailMask) << segmentBits;

    };

    /*
        Specialisation of TwoSegArray.
        User is required to manage de-allocation of memory manually, using the del()
        function, as the deconstructor does not free this space.
        Operations performed on values in the array only modify the "head" segment
        (i.e. the first 32 bits of a double value) unless specified otherwise.
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
            heads = new uint32_t[length];
            tails = new uint32_t[length] ();
        }

        TwoSegArray(uint32_t* heads, uint32_t* tails)
            :heads(heads), tails(tails)
        {}

        ~TwoSegArray() { }

        template<typename T>
        void set(const uint_fast64_t& id, const T& t)
        {
            double d = t;
            doublepun conv;
			conv.d = d;
            heads[id] = (conv.l >> segmentBits);
        }

        template<typename T>
        void setPair(const uint_fast64_t& id, const T& t)
        {
            double d = t;
			doublepun conv;
			conv.d = d;
            tails[id] = (conv.l & tailMask);
            heads[id] = (conv.l >> segmentBits);
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
            heads = new uint32_t[length];
            // we should zero tails when allocating heads array for
            // to avoid unexpected behaviour
            tails = new uint32_t[length] ();
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
        uint32_t* heads;
        uint32_t* tails;

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
        doublepun conv;
		conv.d = d;
        *head = (conv.l >> segmentBits);
        return *this;
    }

    template<typename T>
    inline Head& Head::operator=(const T&& other) noexcept
    {
        double d = other;
        doublepun conv;
		conv.d = d;
        *head = (conv.l >> segmentBits);
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
	   *tail = 0;
        return *this;
    }
    
    template<>
    inline Pair& Pair::operator=(const Head&& other) noexcept
    {
		*head = *other.head;
		*tail = 0;
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
    inline Pair& Pair::operator=(const T&& other) noexcept
    {
        double d = other;
        doublepun conv;
		conv.d = d;
		*tail = (conv.l & tailMask);
		*head = (conv.l >> segmentBits);
        return *this;
    }

    template<typename T>
    inline Pair& Pair::operator=(const T& other)
    {
        double d = other;
        doublepun conv;
		conv.d = d;
		*tail = (conv.l & tailMask);
		*head = (conv.l >> segmentBits);
        return *this;
    }

    template <typename T>
    double Pair::operator+=(const T& rhs)
    {
        double t = *this;
        t += rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator+(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d += rhs;
        return d;
    }

    template<typename T>
    double Pair::operator-=(const T& rhs)
    {
        double t = *this;
        t -= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator-(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d -= rhs;
        return d;
    }

    template<typename T>
    double Pair::operator*=(const T& rhs)
    {
        double t = *this;
        t *= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator*(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d *= rhs;
        return d;
    }

    template<typename T>
    double Pair::operator/=(const T& rhs)
    {
        double t = *this;
        t /= rhs;
        *this = t;
        return t;
    }

    template<typename T>
    inline double operator/(Pair lhs, const T& rhs)
    {
        double d = lhs;
        d /= rhs;
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
