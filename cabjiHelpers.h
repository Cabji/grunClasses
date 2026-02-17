#pragma once
#include <cmath>
#include <functional>
#include <numeric>
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
	double						roundToStep(double val, double step);

	/**
     * @brief Totals a numeric value across a data set, grouped by a specific key.
     * @tparam T The object type (e.g., GrunObject)
     * @tparam KeyFunc A function/lambda that returns the grouping string
     * @tparam ValFunc A function/lambda that returns the numeric value to sum
     */
	 template<typename T> auto totalByGroup(const std::vector<T>& data, auto keySelector, auto valueSelector) 
	{
		// dev-note: template functions MUST be implmented in the .h header file!
		std::unordered_map<std::string, double> report;

		for (const auto& item : data) {
			// keySelector(item) gets the group name (e.g., "32MPa Concrete")
			// valueSelector(item) gets the amount (e.g., 5.5)
			report[keySelector(item)] += valueSelector(item);
		}

		return report;
	}

	std::string	wildcardsToRegexReady(const std::string& str, const bool& useExactSearch = true);

	// Levenshtein Distance (The "Typo" measure)
	// Returns the number of edits to get from a to b
	int levenshtein(const std::string& s1, const std::string& s2);

    // Token Match Score (The "Component" measure)
	/**
	 * @brief Get the token score object (the "Component" measure)
	 * 
	 * @param query 		(std::string)	- 
	 * @param target 		(std::string)	- 
	 * @return 				(double)		- returns 0.0 to 1.0 based on how many words match
	 */
    double get_token_score(const std::string& query, const std::string& target);

    // Simple Containment (The "Subset" measure)
	/**
	 * @brief Get the containment score object (the "Subset" measure)
	 * 
	 * @param query		(std::string)	- 
	 * @param target	(std::string)	- 
	 * @return			(double)		- returns 1.0 if s1 is inside s2, else 0.0
	 */
    double get_containment_score(const std::string& query, const std::string& target);

	/**
	 * @brief Returns the levenshtein distance between 2 strings as a double type
	 * @param s1		(std::string)	- string 1
	 * @param s2		(std::string)	- string 2
	 * @return			(double)		- the resulting score of the difference beteen s1 and s2 (1.0 means they are an exact match)
	 */
	double get_levenshtein_score(const std::string& s1, const std::string& s2);
};
