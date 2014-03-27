#ifndef TOMATL_ENVELOPE_WALKER
#define TOMATL_ENVELOPE_WALKER

namespace tomatl { namespace dsp {

template <typename T> class SimpleWindowedDft
{
public:
	NaiveFt(int length) : mLength(length), mCounter(0)
	{
		mInput = new T[mLength];
		mOutput = new T[mLength];
		mDoubleInput = new T[mLength * 2];
		mPrecomputedWindowFunction = new T[mLength];

		T a0 = 0.35875;
		T a1 = 0.48829;
		T a2 = 0.14128;
		T a3 = 0.01168;

		for (int i = 0; i < mLength; ++i)
		{
			// Blackman-Harris window
			mPrecomputedWindowFunction[i] = (a0 - a1 * std::cos(2 * PI * i / (mLength - 1)) + a2 * cos(4 * PI * i / (mLength - 1)) - a3 * cos(6 * PI * i / (mLength - 1)));
		}
	}

	~NaiveFt()
	{
		delete[] mInput;
		delete[] mOutput;
		delete[] mDoubleInput;
	}

	T* push(T* value)
	{
		T* res = NULL;

		if (mCounter == mLength)
		{
			calculateFast();
			demultiplexOutput();

			mCounter = 0;

			res = mOutput;
		}

		mInput[mCounter] = *value * mPrecomputedWindowFunction[mCounter];
		mDoubleInput[mCounter * 2] = *value * mPrecomputedWindowFunction[mCounter];
		mDoubleInput[mCounter * 2 + 1] = 0.;

		mCounter++;

		return res;
	}

	int getLength()
	{
		return mLength;
	}
private:

	void calculateFast()
		/*
		FFT routine, (C)1996 S.M.Bernsee. Sign = -1 is FFT, 1 is iFFT (inverse)
		Fills fftBuffer[0...2*fftFrameSize-1] with the Fourier transform of the time domain data in fftBuffer[0...2*fftFrameSize-1]. 
		The FFT array takes and returns the cosine and sine parts in an interleaved manner, ie. fftBuffer[0] = cosPart[0], fftBuffer[1] = sinPart[0], asf. 
		fftFrameSize must be a power of 2. It expects a complex input signal (see footnote 2), 
		ie. when working with 'common' audio signals our input signal has to be passed as {in[0],0.,in[1],0.,in[2],0.,...} asf. 
		In that case, the transform of the frequencies of interest is in fftBuffer[0...fftFrameSize].
		*/
	{
		long fftFrameSize = mLength;
		long sign = -1;
		T wr, wi, arg, *p1, *p2, temp;
		T tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
		T* fftBuffer = mDoubleInput;
		long i, bitm, j, le, le2, k, logN;
		logN = (long)(log(fftFrameSize) / log(2.) + .5);

		for (i = 2; i < 2 * fftFrameSize - 2; i += 2)
		{
			for (bitm = 2, j = 0; bitm < 2 * fftFrameSize; bitm <<= 1)
			{
				if (i & bitm) j++;
				j <<= 1;

			}

			if (i < j) 
			{
				p1 = fftBuffer + i; p2 = fftBuffer + j;
				temp = *p1; *(p1++) = *p2;
				*(p2++) = temp; temp = *p1;
				*p1 = *p2; *p2 = temp;
			}
		}

		for (k = 0, le = 2; k < logN; k++)
		{
			le <<= 1;
			le2 = le >> 1;
			ur = 1.0;
			ui = 0.0;
			arg = PI / (le2 >> 1);
			wr = cos(arg);
			wi = sign*sin(arg);

			for (j = 0; j < le2; j += 2)
			{
				p1r = fftBuffer + j; p1i = p1r + 1;
				p2r = p1r + le2; p2i = p2r + 1;

				for (i = j; i < 2 * fftFrameSize; i += le)
				{
					tr = *p2r * ur - *p2i * ui;
					ti = *p2r * ui + *p2i * ur;
					*p2r = *p1r - tr; *p2i = *p1i - ti;
					*p1r += tr; *p1i += ti;
					p1r += le; p1i += le;
					p2r += le; p2i += le;
				}

				tr = ur*wr - ui*wi;
				ui = ur*wi + ui*wr;
				ur = tr;
			}
		}
	}

	void demultiplexOutput()
	{
		for (int i = 0; i < mLength; i += 2)
		{
			mOutput[i / 2] = mDoubleInput[i];
			mOutput[mLength / 2 + i / 2] = mDoubleInput[i + 1];
		}
	}

	// Straightforward implementation of DFT. It may be useful to understand what this thing does,
	// but do not use it for real computations - it's complexity is O(N*2), so with FFT size 1024 (or worse 2048) it will be
	// terribly slow
	void calculate()
	{
		T* mFftSin = mOutput;
		T* mFftCos = mOutput + mLength / 2;

		T arg = 0.;
		T sign = -1.;

		for (int bin = 0; bin < (mLength / 2); ++bin)
		{
			mFftCos[bin] = (mFftSin[bin] = 0.);
			for (int k = 0; k < mLength; ++k)
			{
				arg = 2.0 * (float)bin * PI * (float)k / mLength;
				mFftSin[bin] += mInput[k] * sign * sin(arg);
				mFftCos[bin] += mInput[k] * cos(arg);
			}
		}
	}

	int mLength;
	int mCounter;
	T* mInput;
	T* mDoubleInput;
	T* mOutput;
	T* mPrecomputedWindowFunction;
};
}}
#endif
