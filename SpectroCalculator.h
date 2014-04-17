#ifndef TOMATL_SPECTRO_CALCULATOR
#define TOMATL_SPECTRO_CALCULATOR

namespace tomatl { namespace dsp {

	struct SpectrumBlock
	{
		SpectrumBlock()
		{
			mLength = 0;
			mData = NULL;
			mIndex = 0;
			mSampleRate = 0;
		}

		SpectrumBlock(size_t size, std::pair<double, double>* data, size_t index, size_t sampleRate)
		{
			mLength = size;
			mData = data;
			mIndex = index;
			mSampleRate = sampleRate;
		}

		size_t mLength;
		size_t mIndex;
		size_t mSampleRate;
		std::pair<double, double>* mData;
	};

	template <typename T> class SpectroCalculator
	{
	public:
		SpectroCalculator(std::pair<double, double> attackRelease, size_t index, size_t fftSize = 1024, size_t channelCount = 2)
		{
			mData = new std::pair<double, double>[fftSize];
			mChannelCount = channelCount;
			mFftSize = fftSize;
			mAttackRelease = attackRelease;
			mIndex = index;

			for (int i = 0; i < channelCount; ++i)
			{
				mDfts.push_back(new SimpleWindowedDft<double>(fftSize));
			}
		}

		~SpectroCalculator()
		{
			delete[] mData;

			for (int i = 0; i < mChannelCount; ++i)
			{
				delete mDfts[i];
			}

			mDfts.clear();
		}

		SpectrumBlock process(T* channels, double sampleRate)
		{
			T* response = NULL;

			for (int i = 0; i < mChannelCount; ++i)
			{
				response = mDfts[i]->push(channels + i);
				
			}

			if (response != NULL)
			{
				for (int bin = 0; bin < (mFftSize / 2.); ++bin)
				{
					T ampl = 0.;

					for (int cn = 0; cn < mChannelCount; ++cn)
					{
						T* ftResult = mDfts[cn]->getOutput();
						T* mFftSin = ftResult;
						T* mFftCos = ftResult + mFftSize / 2;

						// http://www.dsprelated.com/showmessage/69952/1.php or see below
						mFftSin[bin] *= 2;
						mFftCos[bin] *= 2;
						mFftSin[bin] /= mFftSize;
						mFftCos[bin] /= mFftSize;

						T nw = std::sqrt(mFftSin[bin] * mFftSin[bin] + mFftCos[bin] * mFftCos[bin]);

						ampl = std::max(nw, ampl);
					}

					double prev = mData[bin].second;

					EnvelopeWalker::staticProcess(&ampl, &prev, &mAttackRelease.first, &mAttackRelease.second, &sampleRate);

					mData[bin].first = bin;
					mData[bin].second = prev;
				}

				return SpectrumBlock(mFftSize / 2., mData, mIndex, sampleRate);
			}
			else
			{
				return SpectrumBlock();
			}
		}

	private:
		std::vector<SimpleWindowedDft<T>*> mDfts;
		std::pair<double, double>* mData;
		std::pair<double, double> mAttackRelease;
		size_t mChannelCount;
		size_t mFftSize;
		size_t mIndex;
	};

}}

/*

What do you expect? I assume you have an N length sinusoidal of
the form

x[n] = A sin(wn)

and find a spectrum coefficient of magnitude A*N/2 (provided
the frquency w is an integer fraction of the sampling frequency
w = m/N, m < N/2).

First the factor 1/2. The spectrum of a real-valued signal is
conjugate symmetric, meaning one real-valed sinusoidal is
represented as two complex-valued sinusoidals according
to Eulers formula,

A*sin(x) = A*(exp(jx)-exp(-jx))/j2.

These two complex-valued sinusoidals have magnitde
A/2, so if you plot only the range [0,fs/2] you need to
scale by a factor 2 to recover A.

Second, the factor N. The DFT is a set of vector products
between the input signal x and N complex elements of
unit magnitude:

X[k] = sum_n=0^N-1 A*exp(j*2*pi*k*n/N)*exp(-j*2*pi*k*n/N)
= N*A

To recover A, one needs to divide by N.

*/

#endif