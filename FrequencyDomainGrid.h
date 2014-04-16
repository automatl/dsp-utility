#ifndef TOMATL_FREQUENCY_DOMAIN_GRID
#define TOMATL_FREQUENCY_DOMAIN_GRID

namespace tomatl{ namespace dsp{

	// TODO: left and top offsets for coordinates.
	class FrequencyDomainGrid
	{
	private:
		tomatl::dsp::OctaveScale mFreqScale;
		tomatl::dsp::LinearScale mMagnitudeScale; // As we use dB values the scale is linear
		tomatl::dsp::Bound2D<double> mBounds;
		tomatl::dsp::Bound2D<double> mFullBounds;
		double* mFreqCache;
		size_t mSampleRate;
		size_t mBinCount;
		size_t mWidth;
		size_t mHeight;
		
		FrequencyDomainGrid(){}
		FrequencyDomainGrid(const FrequencyDomainGrid&){}

		void prepareFreqCache()
		{
			TOMATL_BRACE_DELETE(mFreqCache);

			if (mSampleRate > 0)
			{
				mFreqCache = new double[mSampleRate];
				memset(mFreqCache, 0x0, sizeof(double)* mSampleRate);
			}
		}

		forcedinline int calculateBinNumberToX(double value)
		{
			return mFreqScale.scale(mWidth, mBounds.X, binNumberToFrequency(value), true);
		}

		void recalcGrid()
		{
			mAmplGrid.clear();
			mFreqGrid.clear();

			const int size = 11;
			double freqs[size] = { 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000 };

			for (int i = 0; i < size; ++i)
			{
				if (TOMATL_IS_IN_BOUNDS_INCLUSIVE(freqs[i], mBounds.X.mLow, mBounds.X.mHigh))
				{
					mFreqGrid.push_back(GridLine(freqToX(freqs[i]), freqs[i], freqToString(freqs[i])));
				}
			}

			int dbStep = std::abs(std::abs(mBounds.Y.mHigh) - std::abs(mBounds.Y.mLow)) * 3 < 72 ? 2 : 6;

			for (int i = mBounds.Y.mLow; i <= mBounds.Y.mHigh; i += dbStep)
			{
				mAmplGrid.push_back(GridLine(dbToY(i), i, dbToString(i)));
			}
		}
	public:

		struct GridLine
		{
			GridLine(int location, double value, std::wstring caption)
			{
				mLocation = location;
				mValue = value;
				mCaption = caption;
			}

			int mLocation = 0;
			double mValue = 0.;
			std::wstring mCaption = L"";
		};

		FrequencyDomainGrid(Bound2D<double> fullBounds, size_t sampleRate = 0, size_t binCount = 0, size_t width = 0, size_t height = 0)
		{
			mFullBounds = fullBounds;
			mBounds = fullBounds;
			mSampleRate = sampleRate;
			mBinCount = binCount;
			mWidth = width;
			mHeight = height;
			mFreqCache = NULL;

			prepareFreqCache();
		}

		void updateSize(size_t w, size_t h)
		{
			if (w != mWidth || h != mHeight)
			{
				mHeight = h;
				mWidth = w;
				prepareFreqCache();
				recalcGrid();
			}
		}

		void updateSampleRate(size_t sampleRate)
		{
			if (sampleRate != mSampleRate)
			{
				mSampleRate = sampleRate;
				prepareFreqCache();
			}
		}

		void updateBounds(Bound2D<double> bounds)
		{
			if (!mBounds.areEqual(bounds))
			{
				mBounds = bounds;
				recalcGrid();
				prepareFreqCache();
			}
		}

		void updateBinCount(size_t binCount)
		{
			if (binCount != mBinCount)
			{
				mBinCount = binCount;
			}
		}

		bool isFrequencyVisible(const double& freq) { return TOMATL_IS_IN_BOUNDS_INCLUSIVE(freq, mBounds.X.mLow, mBounds.X.mHigh); }

		size_t getWidth() { return mWidth; }
		size_t getHeight() { return mHeight; }

		size_t getFreqLineCount() { return mFreqGrid.size(); }
		size_t getAmplLineCount() { return mAmplGrid.size(); }

		const GridLine& getFreqLine(int i) { return mFreqGrid[i]; }
		const GridLine& getAmplLine(int i) { return mAmplGrid[i]; }

		forcedinline double binNumberToFrequency(const double& value)
		{
			// Bin count * 2 is just a weird way of getting FFT size
			return value * mSampleRate / (mBinCount * 2);
		}

		forcedinline int freqToX(const double& value)
		{
			if (mFreqCache == NULL)
			{
				return mFreqScale.scale(mWidth, mBounds.X, value, true);
			}

			if (mFreqCache[(int)value] == 0)
			{
				mFreqCache[(int)value] = mFreqScale.scale(mWidth, mBounds.X, value, true);
			}

			return mFreqCache[(int)value];
		}

		
		forcedinline int dbToY(const double& value)
		{
			return mHeight - mMagnitudeScale.scale(mHeight, mBounds.Y, value, true) - 1;
		}

		forcedinline double yToDb(const int& y)
		{
			return mMagnitudeScale.unscale(mHeight, mBounds.Y, mHeight - y - 1, true);
		}

		forcedinline double fullScaleYToDb(const int& y)
		{
			return mMagnitudeScale.unscale(mHeight, mFullBounds.Y, mHeight - y - 1, true);
		}

		forcedinline double fullScaleXToFreq(const double& value)
		{
			return mFreqScale.unscale(mWidth, mFullBounds.X, value, true);
		}

		~FrequencyDomainGrid()
		{
			TOMATL_BRACE_DELETE(mFreqCache);
		}

		static std::wstring freqToString(const double& freq)
		{
			wchar_t buffer[50];
			memset(&buffer, 0x0, 50);

			if (freq >= 1000)
			{
				swprintf(buffer, L"%dk", (int)freq / 1000);
			}
			else
			{
				swprintf(buffer, L"%d", (int)freq);
			}

			return std::wstring(buffer);
		}

		static std::wstring dbToString(const double& ampl)
		{
			wchar_t buffer[50];
			memset(&buffer, 0x0, 50);

			swprintf(buffer, L"%+d", (int)ampl);

			return std::wstring(buffer);
		}
	private:
		std::vector<GridLine> mFreqGrid;
		std::vector<GridLine> mAmplGrid;
	};

}}

#endif