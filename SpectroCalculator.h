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

						// TODO: figure out correct scale
						T nw = 4 * std::sqrt(mFftSin[bin] * mFftSin[bin] + mFftCos[bin] * mFftCos[bin]) / (float)mFftSize;

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

#endif