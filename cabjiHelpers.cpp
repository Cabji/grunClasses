#include "cabjiHelpers.h"
#include <algorithm>
#include <iterator>
#include <sstream>

double	cabji::roundToStep(double val, double step)
{
	if (step == 0) step = 1;
	return std::ceil(val / step) * step;
}

std::string	cabji::wildcardsToRegexReady(const std::string& str, const bool& useExactSearch)
{
	// dev-note: '*' is wildcard character and '?' is match a single character
	std::string	regexStr	= "";

	// convert any wildcard chars (* and ?) to regex equiv.s and escape any other regex significant chars
	for (char c : str) {
		if (c == '*') regexStr += ".*";
		else if (c == '?') regexStr += ".";
		else if (std::string(".+^$|()[]{}").find(c) != std::string::npos) {
			regexStr += "\\"; // Escape standard regex special chars
			regexStr += c;
		}
		else regexStr += c;
	}
	if (useExactSearch) regexStr = "^" + regexStr + "$";
	return regexStr;
}

int cabji::levenshtein(const std::string& s1, const std::string& s2)
{
	const size_t m = s1.size();
	const size_t n = s2.size();

	// Handle empty string cases
	if (m == 0) return static_cast<int>(n);
	if (n == 0) return static_cast<int>(m);

	// We only need two rows of the matrix to calculate the distance
	std::vector<int> prevRow(n + 1);
	std::vector<int> currRow(n + 1);

	// Initialize the first row (representing distances from an empty s1)
	std::iota(prevRow.begin(), prevRow.end(), 0);

	for (size_t i = 1; i <= m; ++i) {
		currRow[0] = static_cast<int>(i);

		for (size_t j = 1; j <= n; ++j) {
			// Cost is 0 if characters match, 1 if they don't
			int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;

			// Calculate minimum of (Deletion, Insertion, Substitution)
			currRow[j] = std::min({
				currRow[j - 1] + 1,      // Insertion
				prevRow[j] + 1,          // Deletion
				prevRow[j - 1] + cost    // Substitution
			});
		}
		// Move current row to previous for next iteration
		prevRow = currRow;
	}

	return prevRow[n];
}

double cabji::get_token_score(const std::string& query, const std::string& target)
{
	auto get_tokens = [](std::string s) 
	{
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		std::stringstream ss(s);
		std::string word;
		std::vector<std::string> tokens;
		while (ss >> word) tokens.push_back(word);
		return tokens;
	};

	auto q_tokens = get_tokens(query);
	auto t_tokens = get_tokens(target);
	
	if (q_tokens.empty()) return 0.0;

	int matches = 0;
	for (const auto& q : q_tokens) 
	{
		for (const auto& t : t_tokens) 
		{
			if (q == t) 
			{
				matches++;
				break; 
			}
		}
	}
	return static_cast<double>(matches) / q_tokens.size();
}

double cabji::get_containment_score(const std::string& query, const std::string& target) 
{
	return (target.find(query) != std::string::npos) ? 1.0 : 0.0;
}

double cabji::get_levenshtein_score(const std::string& s1, const std::string& s2) 
{
	int dist = cabji::levenshtein(s1, s2);
	int maxLength = std::max(s1.length(), s2.length());

	if (maxLength == 0) return 1.0;

	// Convert distance to a similarity percentage (0.0 to 1.0)
	return 1.0 - (static_cast<double>(dist) / maxLength);
}
