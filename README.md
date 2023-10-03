# mantissa-segmentation

Repository for the MEng project "A runtime library for mantissa segmentation".

## Notes about current implementation
- performs mantissa segmentation on a double precision number in two segments of 32 bits
	- "Head" segment - upper 32 bits of a double precision number (1 sign, 11 exponent, 20 mantissa)
	- "Tail" segment - lower 32 bits of a double precision number (32 mantissa)

When compiled with -O3:
	- using "heads" (single segment), performance is generally better than that of a double precision implementation
	  (with some caveats)
    - using pairs (both segments), performance is about 1.5 times as slow as double

Has applications in iterative algorithms that can make use of mixed precision and adaptive precision techniques,
or applications where full double precision range of representation is required, but full 15~ decimal place precision is not.
