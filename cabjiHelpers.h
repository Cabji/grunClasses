#include <cmath>
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

#ifndef CLASS_NAME
#define CLASS_NAME "CabjiHelpers"
#endif

namespace cabji 
{

	/**
	 * @brief rounds a number to the next step
	 * @attention If step == 0, it will be set to = 1 so the function will round to the next whole single unit if it is given 0 as the step
	 * @param val	(double) - the value to round
	 * @param step	(double) - the place/value to round to
	 */
	inline	double	roundToStep(double val, double step)
	{
		if (step == 0) step = 1;
		return std::ceil(val / step) * step;
	}


    /**
     * @brief Totals a numeric value across a data set, grouped by a specific key.
     * @tparam T The object type (e.g., GrunObject)
     * @tparam KeyFunc A function/lambda that returns the grouping string
     * @tparam ValFunc A function/lambda that returns the numeric value to sum
     */
    template<typename T>
    auto totalByGroup(const std::vector<T>& data, 
                      auto keySelector, 
                      auto valueSelector) 
    {
        std::unordered_map<std::string, double> report;

        for (const auto& item : data) {
            // keySelector(item) gets the group name (e.g., "32MPa Concrete")
            // valueSelector(item) gets the amount (e.g., 5.5)
            report[keySelector(item)] += valueSelector(item);
        }

        return report;
    }
}
