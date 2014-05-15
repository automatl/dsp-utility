#ifndef TOMATL_WINDOW_FUNCTION
#define TOMATL_WINDOW_FUNCTION

// TODO: perform some refactoring to avoid implementation inheritance as it would possibly simplify things a little

namespace tomatl { namespace dsp {

class IWindowFunction
{
protected:
	virtual void precalculateWindow(void) = 0;
public:
	// After windowing input signal, it obviously becomes "more quiet", so we're need to compensate that sometimes (http://alpha.science.unitn.it/~bassi/Signal/NInotes/an041.pdf)
	virtual double getWindowScalingFactor() = 0;

	enum SymmetryMode
	{
		modeSymetric = 0,
		modePeriodic
	};
};

template <typename T> class WindowFunction : public IWindowFunction
{
private:
	size_t mLength = 0;
protected:
	T* mPrecalculated = NULL;

	WindowFunction(size_t length)
	{
		mPrecalculated = new T[length];
		mLength = length;
		memset(mPrecalculated, 0x0, sizeof(T)* length);

	}
public:

	forcedinline void applyFunction(T* signal, size_t start, size_t length = 1, bool scale = false)
	{
		int end = start + length;

		int sample = 0;

		for (int i = start; i < end; ++i)
		{
			if (i >= 0 && i < mLength)
			{
				signal[sample] *= mPrecalculated[i];

				if (scale) signal[sample] /= getWindowScalingFactor();
			}
			else
			{
				signal[sample] = 0.;
			}

			++sample;
		}
	}

	// TODO: handle negative indices?
	forcedinline T applyPeriodic(T signal, size_t position)
	{
		T temp = signal;

		applyFunction(&temp, position % getLength());

		return temp;
	}

	forcedinline const size_t& getLength() { return mLength; }

	virtual ~WindowFunction()
	{
		TOMATL_BRACE_DELETE(mPrecalculated);
	}
};

template <typename T> class BlackmanHarrisWindow : public WindowFunction<T>
{
private:
	BlackmanHarrisWindow(const BlackmanHarrisWindow& other) {}
	BlackmanHarrisWindow(BlackmanHarrisWindow &&) {}
	BlackmanHarrisWindow() {}
protected:
	virtual void precalculateWindow()
	{
		T a0 = 0.35875;
		T a1 = 0.48829;
		T a2 = 0.14128;
		T a3 = 0.01168;

		for (int i = 0; i < getLength(); ++i)
		{
			// Blackman-Harris window
			mPrecalculated[i] =
				a0 -
				a1 * std::cos(2 * TOMATL_PI * i / (getLength() - 1)) +
				a2 * std::cos(4 * TOMATL_PI * i / (getLength() - 1)) -
				a3 * std::cos(6 * TOMATL_PI * i / (getLength() - 1));
		}

		int zhopa = 1;
	}
public:
	virtual double getWindowScalingFactor() { return 0.42; }

	BlackmanHarrisWindow(size_t length) : WindowFunction(length) { precalculateWindow(); }
};

template <typename T> class BarlettWindow : public WindowFunction<T>
{
private:
	BarlettWindow(const BarlettWindow& other) {}
	BarlettWindow(BarlettWindow &&) {}
	BarlettWindow() {}

	IWindowFunction::SymmetryMode mMode;
protected:
	virtual void precalculateWindow()
	{
		size_t length = getLength();

		if (mMode == IWindowFunction::modePeriodic) ++length;

		double a = ((double)length - 1.) / 2.;
		double halfL = a;

		for (int i = 0; i < getLength(); ++i)
		{
			mPrecalculated[i] = 1. - std::abs((i - a) / halfL);
		}
	}
public:
	virtual double getWindowScalingFactor() { return 1.; } // TODO: correct scaling

	BarlettWindow(size_t length, IWindowFunction::SymmetryMode mode = IWindowFunction::modeSymetric) : mMode(mode), WindowFunction(length) { precalculateWindow(); }
};

template <typename T> class HannWindow : public WindowFunction<T>
{
private:
	HannWindow(const HannWindow& other) {}
	HannWindow(HannWindow &&) {}
	HannWindow() {}

	IWindowFunction::SymmetryMode mMode;
protected:
	virtual void precalculateWindow()
	{
		size_t length = getLength();

		if (mMode == IWindowFunction::modePeriodic) ++length;

		for (int i = 0; i < getLength(); ++i)
		{
			mPrecalculated[i] = 0.5 * (1. - std::cos(2 * TOMATL_PI * i / (length - 1)));
		}
	}
public:
	virtual double getWindowScalingFactor() { return 0.5; } // TODO: correct scaling

	HannWindow(size_t length, IWindowFunction::SymmetryMode mode = IWindowFunction::modeSymetric) : mMode(mode), WindowFunction(length) { precalculateWindow(); }
};

template <typename T> class SquareWindow : public WindowFunction<T>
{
private:
	SquareWindow(const SquareWindow& other) {}
	SquareWindow(SquareWindow &&) {}
	SquareWindow() {}
protected:
	virtual void precalculateWindow(void)
	{
		for (int i = 0; i < getLength(); ++i)
		{
			mPrecalculated[i] = 1.;
		}
	}
public:
	virtual double getWindowScalingFactor() { return 1.; }

	SquareWindow(size_t length) : WindowFunction(length) { precalculateWindow(); }
};

class WindowFunctionFactory
{
private:
	WindowFunctionFactory(){}
public:
	enum FunctionType
	{
		windowRectangle = 0,
		windowBlackmanHarris,
		windowHann,
		windowBarlett
	};

	template <typename T> static WindowFunction<T>* create(FunctionType type, size_t length, IWindowFunction::SymmetryMode mode = IWindowFunction::modeSymetric)
	{
		if (type == windowRectangle)
		{
			return new SquareWindow<T>(length);
		}
		else if (type == windowBlackmanHarris)
		{
			return new BlackmanHarrisWindow<T>(length);
		}
		else if (type == windowHann)
		{
			return new HannWindow<T>(length, mode);
		}
		else if (type == windowBarlett)
		{
			return new BarlettWindow<T>(length, mode);
		}
		else
		{
			return NULL;
		}
	}
};

}}

#endif