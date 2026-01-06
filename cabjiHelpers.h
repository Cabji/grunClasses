#include <cmath>

#ifndef CLASS_NAME
#define CLASS_NAME "CabjiHelpers"
#endif

namespace cabji 
{

	/**
	 * @brief rounds a number to the next step
	 * @attention If step == 0, it will be converted to = 1 so the function will round to the next whole single unit if it is given 0 as the step
	 * @param val	(double) - the value to round
	 * @param step	(double) - the place/value to round to
	 */
	inline	double	roundToStep(double val, double step)
	{
		if (step == 0) step = 1;
		return std::ceil(val / step) * step;
	}
}