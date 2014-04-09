#ifndef TOMATL_GONIO_CALCULATOR
#define TOMATL_GONIO_CALCULATOR

namespace tomatl { namespace dsp {


template<typename T> class GonioCalculator
{
public:
	GonioCalculator(size_t segmentLength = 512, std::pair<double, double> autoAttackRelease = std::pair<double, double>(0.01, 5000)) : mData(NULL), mProcCounter(0)
	{
		setSegmentLength(segmentLength);
		mSqrt2 = std::pow(2., 0.5);
		mEnvelope.setAttackSpeed(autoAttackRelease.first);
		mEnvelope.setReleaseSpeed(autoAttackRelease.second);
		mEnvelope.setSampleRate(48000); // TODO: listen to sample rate changes somehow (universal solution? interface?)
		
	}

	std::pair<T, T>* handlePoint(const std::pair<T, T>& subject)
	{
		std::pair<T, T> point(subject);
		
		Coord<T>::toPolar(point);
		Coord<T>::rotatePolarDegrees(point, -45.);

		// Scale auto-adjusting. We tend to use max space available even if our signal is not normalized to 0dB
		// TODO: implement custom scale mode (high cap between 0..1)
		double m = std::max(0.01, 1. / mEnvelope.process(point.first));
		// TODO: apply limit to scale auto expantion (when silence is passed, for example) (maybe also controllable)
		point.first *= m;

		Coord<T>::toCartesian(point);

		point.first = std::min(1., std::max(-1., point.first));
		point.second = std::min(1., std::max(-1., point.second));

		mData[mProcCounter] = point;

		++mProcCounter;

		if (mProcCounter >= mSegmentLength)
		{
			mProcCounter = 0;

			return mData;
		}
		else
		{
			return NULL;
		}
	}

	std::pair<T, T>* handlePoint(const T& x, const T& y)
	{
		std::pair<T, T> point(x, y);

		return handlePoint(point);
	}

	GonioCalculator& setSegmentLength(size_t segmentLength)
	{
		TOMATL_DELETE(mData);
		mSegmentLength = segmentLength;

		mData = new std::pair<T, T>[mSegmentLength];

		memset(mData, 0, sizeof(std::pair<T, T>) * mSegmentLength);
		mProcCounter = 0;

		return *this;
	}

	size_t getSegmentLength()
	{
		return mSegmentLength;
	}

	double getCurrentScaleValue()
	{
		return mEnvelope.getCurrentValue();
	}

private:
	size_t mSegmentLength;
	std::pair<T, T>* mData;
	unsigned int mProcCounter;
	double mSqrt2;
	tomatl::dsp::EnvelopeWalker mEnvelope;
};

}}

#endif